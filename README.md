# SDRscan: identify sex chromosomes and the sex-determining region from genome sequencing

**SDRscan is a pipeline to identify sex chromosomes, find the sex-determining
region (SDR), tell an XY system from a ZW system, and determine the genetic sex
of samples, from whole-genome resequencing (WGS) of males and females. It
works for any ploidy (diploid, tetraploid, hexaploid, mixed) and any number of
samples.**

Typical questions it answers:

- *Which chromosome is the sex chromosome in my species?*
- *Where is the sex-determining region, and how large is it?*
- *Is sex determination XY (male heterogametic) or ZW (female heterogametic)?*
- *Which SNPs, k-mers or regions are sex-linked markers I can use for sexing?*
- *Is this individual (or seedling, juvenile, sample of unknown sex) male or female?*
- *Are any of my samples mis-sexed or mislabelled?*

It was built for a dioecious, partly tetraploid plant and is not specific to
any organism: plants, animals, fish, insects or anything with genetic sex
determination and a reference genome.

**Contents:**
[What it does](#what-sdrscan-does) ·
[Requirements](#requirements) ·
[How many samples](#how-many-samples-do-i-need) ·
[Setup](#setup-step-by-step) ·
[Test your installation](#test-your-installation-with-simulated-data) ·
[Configure](#configure) ·
[Run](#run) ·
[Outputs](#outputs) ·
[Interpreting results](#interpreting-results) ·
[Pilot: 1 male + 1 female](#small-studies-one-male-and-one-female-pilot) ·
[Sexing samples](#determining-the-sex-of-samples) ·
[FAQ](#faq-identifying-sex-chromosomes) ·
[Troubleshooting](#troubleshooting) ·
[Citing](#citing)

---

## What SDRscan does

In a species with sex chromosomes, the sex-limited chromosome (**Y** in males
of an XY system, **W** in females of a ZW system) carries sequence and variants
that the other sex lacks. Where X and Y (or Z and W) have stopped recombining,
in the **sex-determining region (SDR)** or sex-linked region, this leaves
signals in genome sequencing data. SDRscan measures them in four independent
ways:

| # | Line of evidence | What it measures | Signal of a sex-linked region |
|---|---|---|---|
| 1 | **Read depth & heterozygosity** (SDRscan's ploidy-aware scan and [findZX](https://github.com/hsigeman/findZX)) | coverage and heterozygosity per sex in genomic windows | the heterogametic sex has lower depth (one copy of X/Z) or different heterozygosity |
| 2 | **Sex-specific SNP density** | SNPs whose ALT allele is carried by (almost) all of one sex and (almost) none of the other | clusters of male-specific (XY) or female-specific (ZW) SNPs |
| 3 | **Sex-specific k-mers** | 31-mers present in most of one sex and absent from every individual of the other; reference-free, then mapped back to the genome | clusters of male- or female-specific sequence |
| 4 | **Assembly comparison** | whole-genome and per-chromosome alignments between assemblies (male vs female, or phased haplotypes) | inversions, deletions or rearrangements on the candidate chromosome |

Each candidate region is then tested:

- **strict vs relaxed thresholds:** high-confidence sites fall in the same region
- **genotype inspection:** the peak is made of many sex-biased SNPs, not a few noisy ones
- **PCA:** samples split by sex in the region; mis-sexed samples stand out
- **sex-label permutation test:** the signal beats random male/female labelling (exact test for small designs)
- **chromosome enrichment test:** the chromosome is enriched against random window sets
- **leave-one-out:** the region survives dropping each individual
- **randomised-sex findZX run:** negative control in which the pattern should disappear

**Male-biased signals indicate XY, female-biased signals indicate ZW.**
SDRscan tests both directions, so the data decide.

Two more modes cover the ends of the sample-size range:

- **Pilot with one male and one female:** depth and heterozygosity only,
  labelled as unreplicated.
- **Genetic sexing:** assigns sex to samples of unknown sex from known
  sex-linked loci, and flags samples whose recorded sex disagrees with their
  genome.

### What SDRscan does not do

- It does not assemble genomes or annotate genes.
- It does not *prove* a region is the SDR. It ranks candidates and measures
  how robust each one is. Agreement between independent lines of evidence is
  the strongest result it gives.
- It does not identify environmental or non-genetic sex determination.
- Gene-level synteny (e.g. GENESPACE riparian plots) is not included.

---

## Requirements

### Data

| Input | Required | Notes |
|---|---|---|
| Paired-end Illumina WGS reads | yes | one pair of FASTQ(.gz) files per individual |
| Sex of each individual | yes (`M`, `F`) or `U` | `U` = unknown: excluded from sex comparisons, can be sexed by stage 13 |
| Reference genome (FASTA, not gzipped) | yes | chromosome-level is best; for an XY test use a female (homogametic) assembly, for ZW a male one |
| Ploidy | yes | `PLOIDY` = 2, 3, 4, 6, ...; per sample if populations are mixed |
| Second assembly / phased haplotypes | optional | for assembly comparison and mapping k-mers to the other sex's genome |

**Sequencing depth (rule of thumb):** about 10× per individual is enough for
read-depth signals. Genotype-based analyses need more, and polyploids more
still: distinguishing 0/0/0/1 from 0/0/1/1 requires several reads per allele
copy. Aim for 20× or more for tetraploids. The simulated tests below use 15×
(diploid) and 20× (tetraploid).

### How many samples do I need?

SDRscan accepts **any number of samples**. Thresholds are fractions of each
sex, so the same settings work for 4 or 400 individuals.

| Analysis (stage) | Minimum | Recommended | Why |
|---|---|---|---|
| QC, alignment, variant calling (1–3) | 1 sample | — | alignment is per sample; variants are called jointly |
| Depth & heterozygosity scan (12) | 1 + 1 (pilot) | ≥ 3 per sex | with ≥ 2 per sex every male must differ from every female |
| findZX (4) | 1 + 1 (pilot) | ≥ 5 per sex | sex means from 1 individual are unreliable |
| Sex-specific SNP density (5) | 2 + 2 | ≥ 5 per sex | with few samples, rare alleles shared by chance look sex-specific |
| Permutation test (6) | 2 + 2 | ≥ 4 per sex | the smallest possible p-value is 1 / C(n, n_focal), see below |
| PCA (6) | 2 + 2 | ≥ 6 in total | descriptive only with fewer |
| Leave-one-out (7) | 3 of the focal sex | ≥ 5 | dropping one must leave at least 2 |
| Sex-specific k-mers (8) | 2 + 2 | ≥ 5 per sex | "present in all but one" needs ≥ 2 of the focal sex |
| Assembly comparison (10) | — | — | uses assemblies, not samples |
| Sexing samples (13) | 1 sample + known sex-linked loci | ≥ 2 known-sex samples per sex | known-sex samples calibrate the scores |

**Smallest achievable p-value of the sex-label permutation test**, whatever
the data:

| Males + females | Distinct labellings | Smallest p |
|---|---|---|
| 2 + 2 | 6 | 0.167 |
| 3 + 3 | 20 | 0.05 |
| 4 + 4 | 70 | 0.014 |
| 5 + 5 | 252 | 0.004 |
| 6 + 6 | 924 | 0.001 |
| 8 + 8 | 12,870 | < 0.001 |

Below about 4 + 4 the test cannot reach p < 0.05. SDRscan still runs it,
reports the floor and warns. When there are few labellings it enumerates them
all (exact test).

Every run starts with a **sample-size report** (counts per sex, unknowns,
ploidy, the p-value floor, available analyses). Stages the design cannot
support are skipped, with the reason given.

### Software

All installed by `environment.yml` with conda:

| Step | Tools |
|---|---|
| reads and alignment | fastp, FastQC, MultiQC, BWA, samtools, GATK4 (read groups only) |
| variant calling and filtering | freebayes, bcftools, htslib (bgzip, tabix) |
| windows and overlaps | bedtools |
| k-mers | KMC3, seqtk |
| assembly comparison | MUMmer4, gnuplot |
| analysis and plots | Python 3 (pandas, numpy, matplotlib, scikit-learn), R (ggplot2, dplyr) |

findZX is included as a git submodule and has its own conda environment.

Operating system: **Linux or macOS. On Windows, use WSL2** (Ubuntu); see
setup step 0.

### Compute

| Dataset | Hardware | Time |
|---|---|---|
| simulated diploid, 2 + 2, 3.3 Mb genome, 15× | 4 threads, 7 GB RAM (laptop, WSL2) | 6 min |
| simulated tetraploid, 6 + 6 + 2 unknown, 3.3 Mb, 20× | 4 threads, 7 GB RAM (laptop, WSL2) | ~31 min |
| real plant genome (~500 Mb), 16 tetraploid samples | 48-core / 250 GB shared server | variant calling: hours to a day per sample; findZX: ~2 weeks on 16 cores |

Memory: KMC uses `KMER_MEM_GB` (default 8 GB; the tests use 3). Other steps
need a few GB; joint-calling memory grows with `MAX_PARALLEL_JOBS` and the
number of samples.
Disk: BAMs are about the size of the FASTQs; the k-mer databases are similar.

---

## Setup (step by step)

**0. Windows only: install WSL2 and Ubuntu.** In PowerShell (as administrator)
run `wsl --install -d Ubuntu-24.04`, restart, open "Ubuntu", and do every
following step inside it. Keep data and the repository in the Linux home
folder (`~`), not under `/mnt/c`, which is much slower.

**1. Install conda** (skip if `conda --version` works). Miniforge is
recommended:

```bash
curl -L -O https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh
bash Miniforge3-$(uname)-$(uname -m).sh -b -p ~/miniforge3
~/miniforge3/bin/conda init bash && exec bash
```

**2. Get SDRscan with its findZX submodule.**

```bash
git clone --recurse-submodules https://github.com/<user>/SDRscan.git
cd SDRscan
# if you cloned without --recurse-submodules:
git submodule update --init
```

**3. Create the environments.**

```bash
conda env create -f environment.yml                             # "sdrscan": everything except findZX
conda env create -f external/findZX/environment.yml -n findZX   # only needed for stage 4
conda activate sdrscan
```

**4. Check the installation** with the simulated test (next section). It
takes a few minutes and runs every stage except findZX.

---

## Test your installation with simulated data

`tests/` simulates a small genome with a known XY sex-determining region,
runs the whole pipeline on it and checks every result against the truth:

```bash
conda activate sdrscan
tests/run_test.sh diploid_2v2    2 2 2          # ploidy 2, 2 males, 2 females
tests/run_test.sh tetraploid_6v6 4 6 6 2 20     # ploidy 4, 6 + 6, 2 of unknown sex, 20x
tests/run_test.sh pilot_1v1      2 1 1          # 1 + 1 pilot
RUN_FINDZX=1 tests/run_test.sh with_findzx 2 3 3   # also stage 4 (slow: builds findZX's environments)
```

The simulated genome (3.3 Mb) has two autosomes, an XY chromosome pair with
a 400 kb SDR, and an unplaced scaffold. The Y carries fixed Y-specific SNPs,
a Y-only insertion, a deletion (lower male depth) and an inversion (visible
in the dotplot). Each test ends with a PASS/FAIL list such as:

```
PASS  stage 5: male-specific SNPs fall in the SDR  (99% in SDR)
PASS  stage 6: permutation test significant  (p=0.00758, floor 0.00108)
PASS  stage 7: SDR survives every leave-one-out  (6/6)
PASS  stage 12: lower male depth over the Y deletion  (1 depth candidates; chr_03:600000-700000 (-0.41))
PASS  stage 13: unknown samples sexed correctly  (U01: F (true F), U02: M (true M))

33/33 checks passed
```

Results of the release tests: diploid 2 + 2 **30/30**, tetraploid 6 + 6 + 2
unknown **33/33**, pilot 1 + 1 **14/14**. On a machine with little memory,
lower the load with `THREADS=4 JOBS=2 KMER_MEM=3 tests/run_test.sh ...`.

Results of the tests run before release are in
[docs/VALIDATION.md](docs/VALIDATION.md). Use `tests/simulate_dataset.py` on
its own to make a dataset with your own number of samples, ploidy and
coverage.

---

## Configure

```bash
cp config/pipeline.env.example config/pipeline.env   # all settings (plain bash, commented)
cp config/samples.example.tsv  config/samples.tsv    # your samples
cp config/pairs.example.tsv    config/pairs.tsv      # assemblies to compare (optional)
```

### Sample sheet: `config/samples.tsv`

Tab-separated, with a header:

```
sample	sex	fq1	fq2	ploidy
ind01	M	/data/ind01_R1.fastq.gz	/data/ind01_R2.fastq.gz	NA
ind02	F	/data/ind02_R1.fastq.gz	/data/ind02_R2.fastq.gz	NA
ind03	U	/data/ind03_R1.fastq.gz	/data/ind03_R2.fastq.gz	NA
```

| Column | Content |
|---|---|
| `sample` | sample ID; becomes the read-group and VCF sample name |
| `sex` | `M`, `F`, or `U` (unknown) |
| `fq1`, `fq2` | absolute paths to the paired read files |
| `ploidy` | optional; overrides `PLOIDY` for this sample (`NA` or empty = use `PLOIDY`) |

### Settings: `config/pipeline.env`

Every setting is commented in the file. The main ones:

| Setting | Meaning | Default |
|---|---|---|
| `REF` | reference FASTA | — |
| `OUTDIR` | output folder | `results` |
| `PLOIDY` | **required**: 2, 3, 4, 6, ... | none |
| `THREADS`, `SORT_MEM` | threads and sort memory per job | 16, 4G |
| `MAX_PARALLEL_JOBS`, `CALL_CHUNK_BP` | joint calling: freebayes jobs at a time, genome chunk size (bp) | 4, 1000000 |
| `MARK_DUPLICATES` | 1 = GATK MarkDuplicatesSpark | 0 |
| `MIN_MAPQ`, `MIN_BASEQ`, `MIN_ALT_COUNT` | freebayes filters | 20, 20, 4 |
| `MIN_QUAL`, `MIN_DP` | merged VCF filters | 30, 5 |
| `CHROM_REGEX` | which sequences are chromosomes (others are ignored in windows and plots) | `chr_?[0-9]+$` |
| `ALLOW_SINGLE_PER_SEX` | 1 = allow a 1 male + 1 female pilot | 0 |
| `FINDZX_RUNS` | findZX runs `SYSTEM\|reference` with SYSTEM = `XY`, `ZW` or `SHUFFLED` | XY, ZW, shuffled |
| `FINDZX_CONDA_ENV`, `FINDZX_BACKGROUND`, `FINDZX_CORES` | findZX environment and resources | `findZX`, 1, 16 |
| `WINDOW`, `STEP` | SNP-density windows (bp) | 100000, 50000 |
| `SNP_RUNS` | named threshold sets `name\|classifier options` | relaxed, strict |
| `RELAXED_RUN`, `STRICT_RUN` | which runs find and validate candidates | relaxed, strict |
| `SNP_DOSAGE_CLASSES`, `SNP_DOSAGE_BASE` | extra runs per ALT dosage (`auto` = 1..ploidy−1) | auto, relaxed |
| `CANDIDATES` | regions to follow up, `chrom:start-end` (fill in after stage 5) | empty |
| `FOCAL_SEX` | `M` to test XY, `F` to test ZW | M |
| `PERMUTATIONS`, `STRONG_ALT`, `STRONG_MAX_OPPOSITE` | permutation test settings (fractions of each sex) | 10000, 0.75, 0.2 |
| `K`, `KMER_MIN_COUNT`, `KMER_MEM_GB` | k-mer length, minimum count, KMC memory (GB; keep below your RAM) | 31, 2, 8 |
| `KMER_RUNS` | `MSK\|n` / `FSK\|n`: in how many samples of the sex a k-mer must occur; n = count, fraction, `all`, `all-1` or `auto` (`all-1` with ≥ 4 samples of that sex, otherwise `all`) | `auto` |
| `KMER_REFS` | `TAG=fasta` for every assembly/haplotype to map k-mers to | — |
| `KMER_WINDOW`, `KMER_STEP`, `KMER_SUBSAMPLE` | k-mer density windows and subsample size | 25000, 5000, 200000 |
| `DEPTH_WINDOW`, `DEPTH_LOG2_THRESHOLD`, `HET_DIFF_THRESHOLD` | depth/heterozygosity scan; `auto` threshold = half the depth drop expected for one missing copy at your ploidy | 100000, auto, 0.1 |
| `SEX_SYSTEM`, `SEXING_REGIONS`, `SEXING_MARKERS`, `SEXING_MSK_DB`, `SEXING_FSK_DB` | sexing samples (stage 13) | empty = use this run's results |
| `DOTPLOT_PAIRS`, `DOTPLOT_MIN_IDENTITY`, `DOTPLOT_MIN_LENGTH`, `DOTPLOT_CHROMS` | assembly comparison | `config/pairs.tsv`, 95, 1000 |

**Testing XY and ZW.** Analyse against a female reference and look for
male-specific signal (XY), and against a male reference and look for
female-specific signal (ZW). Keep two config files with different
`REF`/`OUTDIR` and choose one with `CONFIG=config/other.env ./run_pipeline.sh`.

### Ploidy

`PLOIDY` is required and carried through the whole pipeline:

| Stage | How ploidy is used |
|---|---|
| 2 | joint `freebayes -p <ploidy>` (with `--cnv-map` for mixed ploidy); samples and ploidies are saved in `called_with.tsv` and everything is re-called if they change |
| 3 | refused if the configured ploidy differs from the one used for calling; `sample_ploidy.tsv` is written for later stages |
| 5–7 | genotypes are read at each sample's ploidy; genotypes with the wrong number of alleles are reported and treated as missing; dosage options are checked against the ploidy; dosage-class runs default to 1..ploidy−1 |
| 6 | for ploidy > 2, a second PCA uses ALT dosage instead of presence/absence |
| 12–13 | the depth threshold and the expected depth of the heterogametic sex ((P−1)/P) follow the ploidy |
| 4 | findZX calls variants itself as diploid; for ploidy > 2 its heterozygosity is approximate (depth is unaffected), and SDRscan says so |
| 8 | k-mers do not depend on ploidy |

### SNP thresholds

Each `SNP_RUNS` entry is `name|classifier options`. Add as many as you like:

```bash
SNP_RUNS=(
  "relaxed|--min-called-males 0.75 --min-called-females 0.75 --min-alt-males 0.75 --min-alt-females 0.75 --max-alt-opposite 0.2"
  "strict|--min-called-males 0.85 --min-called-females 0.85 --min-alt-males 0.85 --min-alt-females 0.85 --max-alt-opposite 0"
)
```

| Option | Meaning |
|---|---|
| `--min-called-males`, `--min-called-females` | samples of each sex with a genotype at the site |
| `--min-alt-males`, `--min-alt-females` | samples of the focal sex carrying the ALT allele |
| `--max-alt-opposite` | samples of the other sex allowed to carry it |
| `--alt-dosage N [N ...]` | only genotypes with exactly N ALT copies count (e.g. `1` = 0/0/0/1 in a tetraploid) |
| `--min-alt-dosage N` | genotypes with at least N ALT copies count (default 1) |

Values are counts (`7`) or fractions of that sex (write fractions with a
decimal point). Fractions are rounded towards the stricter side: minimums up
(`0.8` of 6 = 5), maximums such as `--max-alt-opposite` down (`0.25` of 2 =
0). Each run logs the counts the fractions
resolved to, and the number of sites found, in `classifier.log`. Aim for a few
hundred sex-specific sites: a handful means too strict, tens of thousands
means noise. To try one set by hand: `scripts/05_sex_specific_snp_windows.sh
NAME [options]`.

---

## Run

```bash
conda activate sdrscan
./run_pipeline.sh            # stages 1-5 and 12
#  -> look at results/05_snp_windows/*/sex_specific_snp_windows.png and
#     merged_candidate_intervals.tsv, then put the regions in CANDIDATES
./run_pipeline.sh 6 7        # test the candidate regions
./run_pipeline.sh 8 9 10 11  # k-mers, k-mer overlaps, assembly dotplots, figures
./run_pipeline.sh 13         # sex unknown samples, check recorded sexes
./run_pipeline.sh all        # everything (6-7 once CANDIDATES is set; 13 if any sex is U)
CONFIG=config/male_ref.env ./run_pipeline.sh 1 2 3 5   # a second reference
```

Stages skip finished work, so after a failure fix the cause and run the same
command again. Run long jobs in `tmux` or `screen`, or with `nohup`, so they
survive logging out.

| Stage | Script | What it does |
|---|---|---|
| 1 | `scripts/01_qc_align.sh` | fastp trimming, FastQC, `bwa mem`, sort, read groups, MultiQC |
| 2 | `scripts/02_call_variants.sh` | joint freebayes calling of all samples at the declared ploidy, in genome chunks run in parallel |
| 3 | `scripts/03_merge_vcfs.sh` | normalise, quality/depth filter, biallelic SNPs, genotype table |
| 4 | `scripts/04_findzx.sh` | findZX: XY, ZW and shuffled-sex control |
| 5 | `scripts/05_sex_specific_snp_windows.sh` | sex-specific SNPs → windows → ranked candidate intervals |
| 6 | `scripts/06_candidate_followup.sh` | strict overlap, zoom plot, genotype inspection, PCA, permutation tests |
| 7 | `scripts/07_leave_one_out.sh` | does each candidate survive dropping each sample? |
| 8 | `scripts/08_kmer_sex_specific.sh` | male-/female-specific k-mers mapped to every reference |
| 9 | `scripts/09_compare_kmer_runs.sh` | overlap of k-mer candidates between runs |
| 10 | `scripts/10_nucmer_dotplots.sh` | nucmer alignments between assemblies |
| 11 | `scripts/11_figures.sh` | Manhattan plots, dotplots, chromosome correspondence |
| 12 | `scripts/12_depth_het_scan.sh` | ploidy-aware read-depth and heterozygosity scan (also the pilot analysis) |
| 13 | `scripts/13_assign_sex.sh` | assign sex to unknown samples, flag mislabelled ones |

---

## Outputs

Everything goes under `OUTDIR` (default `results/`):

| Folder | Key files |
|---|---|
| `01_align/` | `bam/<sample>.RG.bam`, `fastp/`, `fastqc_*/`, `multiqc/multiqc_report.html` |
| `02_vcf/` | `joint.vcf.gz` (all samples), `called_with.tsv` (sample, ploidy), `chunks/`, `logs/` |
| `03_merged/` | `merged.filtered.snps.vcf.gz`, `merged.genotypes.tsv`, `vcf_samples.txt`, `males.txt`, `females.txt`, `sample_ploidy.tsv` |
| `04_findzx/<run>/` | `units.tsv`, `config.yml`, `findzx.log`, `results/output/no_synteny/{plots,tables}` |
| `05_snp_windows/<run>/` | `sex_specific_sites.tsv`, `*_specific_sites.bed`, `*_specific_windows_<W>_<S>.bed`, `top_candidate_windows.tsv`, **`merged_candidate_intervals.tsv`**, `*_chr_counts.tsv`, **`sex_specific_snp_windows.png`**, `classifier.log` |
| `06_candidates/` | `strict_overlap_counts.tsv`; per region: `*_density.png`, `inspect.*`, `pca_*.png`, **`permutation.result.txt`**, `chrom_enrichment.*` |
| `07_leave_one_out_<region>/` | **`summary.tsv`** (is the region found when each sample is dropped?) |
| `08_kmer/{MSK,FSK}/` | `windows/genome.*_density.<ref>.bed`, `kmer_*_candidates.<ref>.bed`, `*.chrom_counts.tsv` |
| `09_kmer_overlap/` | `overlap_summary.tsv` |
| `10_dotplots/<pair>/` | `*.coords`, `*.png` (mummerplot), `summary.txt` |
| `11_figures/` | `snp_<run>_manhattan.png`, `kmer_*_manhattan.png`, `<pair>_genome.png`, **`<pair>_genome.correspondence_best.tsv`**, per-chromosome dotplots |
| `12_depth_het/` | **`depth_het.candidates.tsv`**, `depth_het.windows.tsv`, `depth_het.genome.png`, `depth_het.per_chrom.tsv` |
| `13_sexing/` | **`sex_assignment.tsv`**, `sex_assignment.png`, `kmer_fractions.tsv` |

The files to look at first are in bold.

---

## Interpreting results

| Result | Strengthens a candidate | Weakens it |
|---|---|---|
| SNP density (5) | several adjacent windows; strict sites inside | isolated windows with 2–3 sites |
| Genotype inspection (6) | many sites with the same sex pattern | the peak rests on 1–2 sites |
| PCA (6) | clean split by sex | no structure; many samples in the "wrong" group |
| Permutation (6) | p < 0.05 and close to the design's floor | p ≥ 0.05 |
| Leave-one-out (7) | region found in every run | disappears when one sample is dropped |
| k-mers (8) | peaks on the same chromosome or region | peaks scattered genome-wide |
| Depth scan (12), findZX (4) | lower depth / different heterozygosity in one sex over several adjacent windows (`support` = `strong`); gone in the shuffled run | single `weak` windows (expected by chance, especially in pilots); same pattern in the shuffled run |
| Dotplots (10–11) | structural change on the candidate chromosome | absence is not evidence against |

**XY or ZW?** Male-specific SNPs and k-mers with a female reference, and lower
male depth over X-like regions, indicate XY. Female-specific signals with a
male reference, and lower female depth, indicate ZW.

Check the chromosome correspondence table (stage 11) before comparing results
between assemblies: chromosome numbers do not always refer to the same
chromosome in different assemblies.

---

## Small studies: one male and one female (pilot)

With one individual of a sex there is no replication. A difference could be
between those two individuals rather than between the sexes, and no test is
possible (2 labellings, so p ≥ 0.5). SDRscan refuses sex comparisons unless
`ALLOW_SINGLE_PER_SEX=1`, and then runs only analyses that still mean
something:

| Runs in a pilot | Why |
|---|---|
| depth scan (12), findZX (4) | a differentiated X/Z region gives lower depth in the heterogametic sex, and sex-limited sequence a larger drop; one individual shows this |
| heterozygosity scan (12) | hemizygous regions lose heterozygosity (noisy with one individual) |
| assembly comparison (10) | uses assemblies, not samples |

| Skipped in a pilot | Why |
|---|---|
| SNP density, k-mers (5, 8, 9) | with n = 1, every private variant looks "sex-specific" |
| permutation, PCA, leave-one-out (6, 7) | need replication |

Pilot outputs are labelled `PILOT`. Use them to decide where to look and how
many samples to sequence. Young sex chromosomes with a small SDR are usually
invisible in a pilot.

---

## Determining the sex of samples

Mark samples of unknown sex as `U` in the sample sheet and run stage 13.
Each sample is scored on up to three kinds of evidence at known sex-linked
loci:

| Evidence | Uses | Male-like when (XY example) |
|---|---|---|
| depth | depth over `SEXING_REGIONS` (default: the stage-12 depth candidates, else `CANDIDATES`) | lower depth over X-like regions ((P−1)/P of normal) |
| snp | genotypes at `SEXING_MARKERS` (default: the strict stage-5 sites) | carries male-specific ALT alleles |
| kmer | share of the stage-8 male- and female-specific k-mer sets in the sample | contains male-specific k-mers |

- **Calibration:** scores are calibrated on the samples of known sex (typical
  male ≈ +1, typical female ≈ −1). With no known-sex samples, theoretical
  expectations are used, and depth then needs `SEX_SYSTEM=XY` or `ZW`.
- **Calls:** M or F with high, medium or low confidence, or *uncertain* when
  the evidence is weak or conflicting.
- **Mislabel check:** samples of recorded sex are re-checked, and
  disagreements are flagged `MISMATCH`.
- **Circularity:** markers discovered from the same samples separate them
  trivially, so for known-sex samples depth is the independent check.
- **Markers from an earlier project:** set `SEXING_REGIONS`,
  `SEXING_MARKERS`, `SEXING_MSK_DB` and `SEXING_FSK_DB` to that project's
  results to sex new samples of the same species without re-running
  discovery.
- **No known locus:** sexing is not possible. Two unsexed individuals cannot
  be told apart from genome data alone.

---

## FAQ: identifying sex chromosomes

**How do I identify sex chromosomes from genome sequencing data?**
Sequence several males and females (whole-genome resequencing, ideally ≥ 5
per sex), map them to a reference genome and compare the sexes. Sex
chromosomes show (1) lower read depth in the heterogametic sex over
hemizygous regions, (2) SNPs present in one sex only, and (3) sequence
(k-mers) present in one sex only. SDRscan runs all three and tests the
candidate regions.

**How do I find the sex-determining region (SDR)?**
Run stages 1–5 and 12, and look for clusters of sex-specific SNPs and depth
differences. Then confirm the candidate with stages 6–8: strict sites,
permutation test, leave-one-out and k-mers.

**How can I tell whether a species is XY or ZW?**
In an XY system, males carry the Y-specific variants: male-specific SNPs and
k-mers, and lower male depth over X-only regions. In a ZW system the same
signals appear in females. SDRscan reports both directions.

**Can I identify sex chromosomes without a reference genome?**
The k-mer step (stage 8) finds sex-specific sequence without a reference. A
reference is needed to place it on chromosomes and for the other analyses.

**How many males and females do I need to sequence?**
At least 2 per sex to compare the sexes, and 4–5 or more per sex for
statistical tests that can reach p < 0.05. See [How many samples do I
need?](#how-many-samples-do-i-need).

**Does it work for polyploids?**
Yes. Ploidy is a required setting and is used at every step. It was developed
on a tetraploid plant and tested on simulated diploid and tetraploid data.

**How do I find sex-linked markers for sexing plants or animals?**
The strict sex-specific SNPs (`05_snp_windows/strict/sex_specific_sites.tsv`)
and the sex-specific k-mers (stage 8) are candidate sex-linked markers. Stage
13 uses them to sex new samples.

**How do I determine the sex of a sample (genetic or molecular sexing)?**
Once sex-linked loci are known, mark the sample `U` and run stage 13.

**I only have one male and one female. What can I do?**
Run a pilot (`ALLOW_SINGLE_PER_SEX=1`); see [Small
studies](#small-studies-one-male-and-one-female-pilot).

**Why are variants called jointly rather than one sample at a time?**
A single-sample VCF only lists sites where that sample differs from the
reference. After merging, a sample that matches the reference at a site
appears as *missing*, not as 0/0. At a Y-specific SNP every female then looks
missing, and the site fails the "called in enough females" filter: the
sex-specific sites disappear. In the simulated test, per-sample calling found
0 male-specific SNPs; joint calling recovers them. If you have per-sample
VCFs from elsewhere, re-call jointly (or use gVCFs) before stage 3.

**How is SDRscan different from findZX?**
findZX (included) uses read depth and heterozygosity, assuming diploidy.
SDRscan adds ploidy-aware genotype analyses, sex-specific SNP density,
reference-free k-mers, assembly comparison, statistical validation of
candidates, small-sample handling and sexing of unknown samples.

---

## Troubleshooting

| Problem | Fix |
|---|---|
| `PLOIDY must be set` | set `PLOIDY` in `config/pipeline.env` |
| `sex must be M, F or U` | the sample sheet's `sex` column has another value (e.g. `male`) |
| `need at least 2 males and 2 females` | add samples, or set `ALLOW_SINGLE_PER_SEX=1` for a pilot |
| `samples or ploidy differ from those used for calling` | you changed `PLOIDY`, a `ploidy` column or the sample list: re-run stage 2 (it re-calls everything) |
| `WARNING: N genotypes did not match the declared ploidy` | the VCF was called at another ploidy, or a sample's `ploidy` column is wrong |
| `only N of M chunks finished` | see `02_vcf/logs/<chunk>.log`; often a truncated BAM (delete it and re-run stage 1) or too little memory (lower `MAX_PARALLEL_JOBS`). Re-running stage 2 resumes from the finished chunks |
| `samples missing from the VCF` | read-group sample names must equal the sample IDs; re-run stage 1 for those samples |
| `fastqc: unknown option -o` or similar | an old copy of SDRscan; FastQC ≥ 0.13 needs `--outdir`, which current SDRscan uses |
| `conda env create` takes very long | use Miniforge (conda-forge only) rather than Anaconda's `defaults` channel, or `mamba env create -f environment.yml` |
| KMC fails or is killed | lower `KMER_MEM_GB` below the free RAM, or use a disk with more space for `tmp/` |
| far too many or too few sex-specific SNPs | change `SNP_RUNS` thresholds; check `classifier.log` |
| findZX fails while building environments | keep `FINDZX_CONDA_PREFIX` on a local disk and re-run; findZX resumes |
| findZX plots show fewer samples than expected | check `04_findzx/<run>/units.tsv`; it is rebuilt from the sample sheet each run |
| a shared server becomes slow or unresponsive | lower `THREADS`, `MAX_PARALLEL_JOBS` and `FINDZX_CORES` |
| the job stops when you log out | run inside `tmux`/`screen` or with `nohup` |
| bedtools `WARNING: ... inconsistent naming convention` | harmless: the genome mixes names such as `chr_01` and `scaffold_1` |

---

## Standalone tools

Every analysis step is also a script you can use on its own (`--help` for
options):

| Tool | Purpose |
|---|---|
| `bin/classify_sex_specific_snps.py` | male-/female-specific SNPs from a genotype table, any ploidy |
| `bin/plot_sex_specific_windows.py` | genome-wide window plot + top windows |
| `bin/merge_top_candidate_windows.py` | merge top windows into ranked candidate intervals |
| `bin/strict_overlap_with_candidates.py` | strict sites inside relaxed candidates |
| `bin/inspect_candidate_region.py`, `bin/run_candidate_region_inspection.py` | per-site sex bias in a region / top-N regions |
| `bin/plot_candidate_chromosome.py` | zoomed plot of one chromosome |
| `bin/pca_from_gt_presence.py` | PCA of a region (presence or dosage) |
| `bin/permute_sex_labels.py` | sex-label permutation test (exact when small) |
| `bin/window_enrichment_permutation.py` | chromosome enrichment against random window sets |
| `bin/depth_het_scan.py` | read depth and heterozygosity between sexes |
| `bin/assign_sex.py` | genetic sexing from depth, SNP and k-mer evidence |
| `bin/plot_manhattan_windows.py` | several runs/references on one Manhattan plot |
| `findzx/make_units.py` | findZX sample sheet for XY, ZW or shuffled sexes |
| `tests/simulate_dataset.py` | simulate a dataset with a known SDR |

## Repository layout

```
run_pipeline.sh     runs the stages in order (prints a sample-size report first)
scripts/            stages 01-13 and common.sh (config, sample sheet, ploidy, sample-size checks)
bin/                Python analysis and plotting tools (gtutils.py = genotype parser)
R/                  dotplots and chromosome correspondence
findzx/             findZX helpers
external/findZX/    findZX (git submodule, unmodified)
config/             example settings, sample sheet, assembly pairs
tests/              simulator, end-to-end test, result checks
docs/WORKFLOW.md    how the analysis was developed, step by step
docs/VALIDATION.md  results of the simulated end-to-end tests
```

## Citing

If you use SDRscan, please cite it (see `CITATION.cff`, or "Cite this
repository" on GitHub) and the tools it runs: findZX (see
`external/findZX/README.md`), BWA, freebayes, bcftools/htslib, bedtools, KMC,
MUMmer, fastp and samtools.

**Keywords:** sex chromosome identification, how to identify sex
chromosomes, sex-determining region, SDR, sex-linked region, XY vs ZW, sex
determination system, sex-specific SNPs, sex-specific k-mers, sex-linked
markers, genetic sexing, molecular sexing, dioecious plants, polyploid
genomics, whole-genome resequencing, read depth, heterozygosity, findZX,
bioinformatics pipeline.
