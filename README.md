# SDRscan

**Find the sex-determining region of a genome from resequenced males and females, and tell XY from ZW.**

SDRscan takes whole-genome sequencing reads from individuals of known sex and
a reference genome, and reports which chromosome regions differ consistently
between the sexes. It combines four independent lines of evidence, then tests
each candidate region before you trust it. It works for any ploidy and any
number of samples from 2 males + 2 females upwards.

---

## What SDRscan does

In a sex-chromosome system the sex-limited chromosome (Y in males, W in
females) carries sequence and variants that the other sex does not have.
SDRscan looks for that signal in four ways:

| # | Line of evidence | What it measures | Signal of a sex-linked region |
|---|---|---|---|
| 1 | **Depth & heterozygosity** (via [findZX](https://github.com/hsigeman/findZX)) | read depth and heterozygosity per sex in windows | depth or heterozygosity differs between sexes |
| 2 | **Sex-specific SNP density** | SNPs whose ALT allele is carried by (almost) all of one sex and (almost) none of the other | windows enriched for such SNPs |
| 3 | **Sex-specific k-mers** | 31-mers present in most of one sex and absent from every individual of the other, mapped back to the genome | windows enriched for such k-mers |
| 4 | **Assembly comparison** | whole-genome / per-chromosome alignments between assemblies (e.g. male vs female, or two haplotypes) | structural differences on the candidate chromosome |

Candidate regions from (2) are then checked:

- **strict vs relaxed thresholds:** strict sites fall inside relaxed candidates
- **genotype inspection:** many repeated sex-biased SNPs, not a few noisy ones
- **PCA:** samples split by sex, and mis-sexed samples show up
- **sex-label permutation test:** is the signal stronger than under random labels?
- **chromosome enrichment test:** is the candidate chromosome enriched against random window sets?
- **leave-one-out:** does the candidate survive dropping each sample?
- **randomised-sex findZX run:** negative control; the pattern should disappear

**Male-biased signals point to XY, female-biased to ZW.** SDRscan runs both
directions so the data decide.

### What it does not do

- It does not assemble genomes or annotate genes in the candidate region.
- It does not prove a region is the SDR. It ranks candidates and quantifies
  how robust each one is. Agreement between independent lines of evidence is
  the strongest result it can give.
- Gene-level synteny between assemblies (e.g. GENESPACE riparian plots) is
  not included.

---

## Requirements

### Data

| Input | Needed | Notes |
|---|---|---|
| Paired-end WGS reads (Illumina) | yes | one pair of FASTQ files per individual |
| Sex of every individual | yes | phenotypic sex; the PCA can flag likely mis-sexed samples |
| Reference genome (FASTA) | yes | chromosome-level is best; scaffold-level works but plots are harder to read |
| Ploidy | yes | set `PLOIDY` (2, 3, 4, 6, ...), per sample if populations are mixed |
| Second assembly / haplotypes | optional | needed for assembly comparison and for mapping k-mers to the other sex's genome |

**Sequencing depth (rule of thumb).** About 10× per individual is enough for
read-depth signals. Genotype-based analyses need more, and polyploids more
still, because distinguishing 0/0/0/1 from 0/0/1/1 requires several reads per
allele copy. Aim for 20× or more for tetraploids.

### How many samples do I need?

SDRscan accepts **any number of samples**. Thresholds are given as fractions
of each sex, so the same settings work for 4 or 400 individuals. Analyses that
compare the sexes require **at least 2 males and 2 females**, and some
analyses need more:

| Analysis | Minimum | Recommended | Why |
|---|---|---|---|
| QC, alignment, variant calling (stages 1–3) | 1 sample | — | per-sample steps |
| findZX (stage 4) | 2 + 2 | ≥ 5 per sex | findZX itself runs with 1 + 1, but sex means from 1 individual are unreliable |
| Sex-specific SNP density (stage 5) | 2 + 2 | ≥ 5 per sex | with few samples, chance sharing of rare alleles looks "sex-specific" |
| Sex-label permutation test (stage 6) | 2 + 2 | ≥ 4 per sex | the smallest possible p-value is 1 / C(n, n_focal), see below |
| PCA (stage 6) | 2 + 2 | ≥ 6 in total | descriptive only with fewer |
| Leave-one-out (stage 7) | 3 of the focal sex | ≥ 5 | dropping one must leave at least 2 |
| Sex-specific k-mers (stage 8) | 2 + 2 | ≥ 5 per sex | "present in all but one" needs ≥ 2 of the focal sex |
| Assembly comparison (stage 10) | — | — | uses assemblies, not samples |

**Smallest achievable p-value of the sex-label permutation test**, whatever
the data:

| Design (males + females) | Distinct labellings | Smallest p |
|---|---|---|
| 2 + 2 | 6 | 0.167 |
| 3 + 3 | 20 | 0.05 |
| 4 + 4 | 70 | 0.014 |
| 5 + 5 | 252 | 0.004 |
| 8 + 8 | 12,870 | < 0.001 |

With fewer than about 4 + 4 the test cannot reach conventional significance.
SDRscan then still runs it, reports the floor, and warns. When the number of
labellings is small it enumerates all of them (exact test); otherwise it
samples `PERMUTATIONS` of them.

Every run starts with a **sample-size report**: counts per sex, ploidy, the
permutation p-value floor, and which analyses your design supports.

### Software

Installed by `environment.yml` (conda): fastp, FastQC, MultiQC, BWA, samtools,
GATK4 (read groups), freebayes, bcftools/htslib, bedtools, KMC3, seqtk,
MUMmer4 + gnuplot, Python 3 (pandas, numpy, matplotlib, scikit-learn), R
(ggplot2, dplyr). findZX comes as a git submodule and uses its own conda
environment. Linux or macOS. On Windows, use WSL.

### Compute

Read alignment, variant calling and findZX dominate run time. As a reference
point, the original analysis ran on a shared 48-core / 250 GB server:

- per-sample variant calling took hours to a day per sample
- one findZX run on 16 samples took about two weeks on 16 cores

KMC needs a few GB of RAM and fast scratch disk. Plan disk space for BAMs
(roughly the size of the FASTQs) and the k-mer databases.

---

## Install

```bash
git clone --recurse-submodules https://github.com/<user>/SDRscan.git
cd SDRscan
conda env create -f environment.yml                               # SDRscan tools
conda env create -f external/findZX/environment.yml -n findZX     # findZX (stage 4)
conda activate sdrscan
```

If you cloned without `--recurse-submodules`: `git submodule update --init`.

## Configure

```bash
cp config/pipeline.env.example config/pipeline.env   # settings (plain bash)
cp config/samples.example.tsv  config/samples.tsv    # your samples
cp config/pairs.example.tsv    config/pairs.tsv      # assemblies to compare (optional)
```

**`config/samples.tsv`** (tab-separated, header required):

| column | content |
|---|---|
| `sample` | sample ID (used as the VCF sample name) |
| `sex` | `M` or `F` |
| `fq1`, `fq2` | paths to the read files |
| `ploidy` | optional; overrides `PLOIDY` for this sample (`NA` = use `PLOIDY`) |

**`config/pipeline.env`**, main settings:

| Setting | Meaning |
|---|---|
| `REF`, `OUTDIR` | reference FASTA and output folder |
| `PLOIDY` | **required**, no default |
| `THREADS`, `MAX_PARALLEL_SAMPLES` | resources (keep modest on shared servers) |
| `FINDZX_RUNS` | findZX runs: `XY|female_ref.fa`, `ZW|male_ref.fa`, `SHUFFLED|...` |
| `SNP_RUNS` | named threshold sets for the SNP analysis (see below) |
| `SNP_DOSAGE_CLASSES` | `auto` = one extra run per ALT dosage 1..ploidy−1 |
| `CANDIDATES`, `FOCAL_SEX` | regions to follow up (filled in after stage 5) and the sex to test |
| `KMER_RUNS`, `KMER_REFS` | k-mer modes (`MSK`/`FSK` with `all-1`, a count or a fraction) and references to map to |
| `CHROM_REGEX` | which sequences count as chromosomes in windows and plots |

To test both XY and ZW properly, analyse against a female reference (look
for male-specific signal) and a male reference (female-specific signal). Keep
two config files with different `REF`/`OUTDIR` and pass one with `CONFIG=...`.

### Ploidy

`PLOIDY` is required and carried through every step: freebayes calls at that
ploidy and records it per sample (stage 2); stage 3 refuses to merge VCFs
called at a different ploidy than configured and writes `sample_ploidy.tsv`;
every genotype-based tool reads genotypes at that ploidy, reports genotypes
with the wrong number of alleles, and checks dosage options against it.

findZX does its own diploid variant calling. For ploidy > 2 its
heterozygosity output is approximate (its read-depth output is unaffected),
and SDRscan prints a note when this applies.

### SNP thresholds

Each `SNP_RUNS` entry is `name|classifier options`. Add as many as you like:

```bash
SNP_RUNS=(
  "relaxed|--min-called-males 0.75 --min-called-females 0.75 --min-alt-males 0.75 --min-alt-females 0.75 --max-alt-opposite 2"
  "strict|--min-called-males 0.85 --min-called-females 0.85 --min-alt-males 0.85 --min-alt-females 0.85 --max-alt-opposite 0"
)
```

| Option | Meaning |
|---|---|
| `--min-called-males/-females` | samples of that sex with a genotype |
| `--min-alt-males/-females` | samples of the focal sex carrying ALT |
| `--max-alt-opposite` | samples of the other sex allowed to carry ALT |
| `--alt-dosage N [N ...]` | only genotypes with exactly N ALT copies count as ALT-present |
| `--min-alt-dosage N` | genotypes with ≥ N ALT copies count (default 1) |

Values are counts (`7`) or fractions of that sex (`0.8` = 80%, rounded up).
Each run logs the counts it resolved to and the number of sites found in
`classifier.log`. Aim for a few hundred sex-specific sites: very few means too
strict, tens of thousands means noise. You can also run one set directly:
`scripts/05_sex_specific_snp_windows.sh NAME [options]`.

---

## Run

```bash
./run_pipeline.sh            # stages 1-5
#   -> look at results/05_snp_windows/*/sex_specific_snp_windows.png and
#      merged_candidate_intervals.tsv, then set CANDIDATES in the config
./run_pipeline.sh 6 7        # test the candidates
./run_pipeline.sh 8 9 10 11  # k-mers, overlaps, assembly dotplots, figures
./run_pipeline.sh all        # everything (6-7 only once CANDIDATES is set)
```

Stages skip finished work, so after a failure fix the cause and run again.
Run long jobs in `tmux` or with `nohup`.

| Stage | Script | What it does | Main output (`results/`) |
|---|---|---|---|
| 1 | `01_qc_align.sh` | fastp, FastQC, bwa mem, sort, read groups, MultiQC | `01_align/bam/` |
| 2 | `02_call_variants.sh` | freebayes per sample at the declared ploidy | `02_vcf/` |
| 3 | `03_merge_vcfs.sh` | normalise, merge, QUAL ≥ 30 & DP ≥ 5, biallelic SNPs, genotype table | `03_merged/` |
| 4 | `04_findzx.sh` | findZX: XY, ZW and shuffled-sex control | `04_findzx/<run>/results/` |
| 5 | `05_sex_specific_snp_windows.sh` | classify sex-specific SNPs, count in windows, rank candidate intervals | `05_snp_windows/<run>/` |
| 6 | `06_candidate_followup.sh` | strict overlap, zoom, genotype inspection, PCA, permutation tests | `06_candidates/<region>/` |
| 7 | `07_leave_one_out.sh` | drop each focal-sex sample; does the region persist? | `07_leave_one_out_<region>/` |
| 8 | `08_kmer_sex_specific.sh` | male-/female-specific k-mers mapped to each reference | `08_kmer/{MSK,FSK}/` |
| 9 | `09_compare_kmer_runs.sh` | overlap of k-mer candidates between runs | `09_kmer_overlap/` |
| 10 | `10_nucmer_dotplots.sh` | nucmer alignments between assemblies | `10_dotplots/` |
| 11 | `11_figures.sh` | Manhattan plots, dotplots, chromosome correspondence table | `11_figures/` |

**[docs/WORKFLOW.md](docs/WORKFLOW.md)** explains each step in the order the
analysis was developed: the question it answers, the settings, what the
results showed, and why the next step was needed.

## Interpreting results

| Result | Strengthens a candidate | Weakens it |
|---|---|---|
| SNP density | several adjacent windows; strict sites inside | single windows with 2–3 sites |
| Genotype inspection | many sites with the same sex pattern | the peak rests on 1–2 sites |
| PCA | clean split by sex | no structure, or many samples in the "wrong" cluster |
| Permutation | p small *and* below the design's floor | p near the floor or above 0.05 |
| Leave-one-out | candidate found in every run | disappears when one sample is dropped |
| k-mers | peaks on the same chromosome/region | peaks scattered genome-wide |
| findZX | depth/heterozygosity shift on the same chromosome; gone in the shuffled run | same pattern in the shuffled run |
| Dotplots | structural change on the candidate chromosome | — (absence is not evidence against) |

Before comparing results across assemblies, check the chromosome
correspondence table from stage 11: chromosome numbers do not always mean the
same chromosome in different assemblies.

## Standalone tools (`bin/`)

Every analysis step is also a standalone script (`--help` for options):
`classify_sex_specific_snps.py`, `plot_sex_specific_windows.py`,
`merge_top_candidate_windows.py`, `strict_overlap_with_candidates.py`,
`inspect_candidate_region.py`, `run_candidate_region_inspection.py`,
`plot_candidate_chromosome.py`, `pca_from_gt_presence.py`,
`permute_sex_labels.py`, `window_enrichment_permutation.py`,
`plot_manhattan_windows.py`. `gtutils.py` is the shared genotype parser.
`findzx/make_units.py` writes findZX sample sheets (XY, ZW or `--shuffle`).

## Layout

```
run_pipeline.sh   runs the stages in order (prints a sample-size report first)
scripts/          stages 01-11 + common.sh (config, sample sheet, ploidy, sample-size checks)
bin/              Python analysis and plotting tools
R/                dotplots and chromosome correspondence
findzx/           findZX helpers
external/findZX   findZX (git submodule, unmodified)
config/           example settings, sample sheet, assembly pairs
docs/WORKFLOW.md  step-by-step account of the analysis and its decisions
```

## Citing

If you use SDRscan, please cite this repository and the tools it runs: findZX
(see `external/findZX/README.md` for its citation), BWA, freebayes,
bcftools/htslib, bedtools, KMC, MUMmer, fastp and samtools.
