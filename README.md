# sex-region-pipeline

A pipeline for locating the sex-determining / sex-linked region of a genome
from whole-genome resequencing of males and females, and for telling XY from
ZW. It combines independent lines of evidence and checks each candidate
before trusting it:

| Evidence | Idea | Stage |
|---|---|---|
| **Depth & heterozygosity** ([findZX](https://github.com/hsigeman/findZX)) | a sex-linked region differs between sexes in read depth and heterozygosity | 4 |
| **Sex-specific SNP density** | SNPs whose ALT allele is (almost) only found in one sex, counted in sliding windows | 5 |
| **Candidate checks** | strict-site overlap, genotype inspection, PCA, permutation tests, leave-one-out | 6–7 |
| **Sex-specific k-mers** | sequence present in most of one sex and absent from all of the other, mapped back to each assembly | 8–9 |
| **Assembly comparison** | whole-genome / per-chromosome dotplots and a chromosome correspondence table between assemblies | 10 |

Genotypes are reduced to ALT-present / ALT-absent per sample, so the SNP
analyses work for diploids and polyploids (it was built on tetraploid calls).

**[docs/WORKFLOW.md](docs/WORKFLOW.md)** walks through the analysis in the
order it was done: why each step was added, what it showed, how thresholds
were chosen, and the lessons learned.

---

## Install

```bash
git clone --recurse-submodules https://github.com/<user>/sex-region-pipeline.git
cd sex-region-pipeline
conda env create -f environment.yml          # tools for stages 1-3, 5-11
conda env create -f external/findZX/environment.yml -n findZX   # stage 4
conda activate sex-region
```

findZX is included as a git submodule (`external/findZX`, pinned commit) and
used unmodified. If you cloned without `--recurse-submodules`, run
`git submodule update --init`.

## Configure

```bash
cp config/pipeline.env.example config/pipeline.env    # paths, threads, ploidy, run lists
cp config/samples.example.tsv  config/samples.tsv     # sample <tab> sex (M/F) <tab> fq1 <tab> fq2
cp config/pairs.example.tsv    config/pairs.tsv       # assemblies to compare (stage 10)
```

`pipeline.env` is plain bash. Besides paths and resources it holds the list
of runs for each stage: `FINDZX_RUNS`, `SNP_RUNS` (threshold sets),
`CANDIDATES`, `KMER_RUNS` and `KMER_REFS`. To analyse against two references
(female assembly for XY, male assembly for ZW), keep two config files with
different `REF`/`OUTDIR` and run each with `CONFIG=...`.

## Run

```bash
./run_pipeline.sh            # stages 1-5: align, call, merge, findZX, SNP windows
# look at results/05_snp_windows/*/sex_specific_snp_windows.png and
# merged_candidate_intervals.tsv, then set CANDIDATES in the config
./run_pipeline.sh 6 7        # follow up the candidates
./run_pipeline.sh 8 9 10 11  # k-mers, k-mer overlaps, dotplots, figures
./run_pipeline.sh all        # everything (6-7 only if CANDIDATES is set)
```

Stages skip outputs that already exist, so a failed run can be resumed. Run
long jobs in `tmux` or with `nohup`; on shared servers keep `THREADS` and
`MAX_PARALLEL_SAMPLES` modest.

| Stage | Script | Does | Main outputs (`results/`) |
|---|---|---|---|
| 1 | `01_qc_align.sh` | fastp, FastQC, bwa mem, sort, read groups, MultiQC | `01_align/bam/*.RG.bam` |
| 2 | `02_call_variants.sh` | freebayes per sample (`PLOIDY`) | `02_vcf/*.vcf.gz` |
| 3 | `03_merge_vcfs.sh` | norm, rename, merge, QUAL/DP filter, biallelic SNPs | `03_merged/merged.genotypes.tsv` |
| 4 | `04_findzx.sh` | findZX XY / ZW / shuffled-sex runs | `04_findzx/<run>/results/` |
| 5 | `05_sex_specific_snp_windows.sh` | classify sites, window counts, candidate intervals, plot | `05_snp_windows/<run>/` |
| 6 | `06_candidate_followup.sh` | strict overlap, zoom, genotype inspection, PCA, permutations | `06_candidates/<region>/` |
| 7 | `07_leave_one_out.sh` | drop each focal-sex sample, does the region persist? | `07_leave_one_out_<region>/summary.tsv` |
| 8 | `08_kmer_sex_specific.sh` | KMC male/female-specific k-mers mapped to each reference | `08_kmer/{MSK,FSK}/windows/` |
| 9 | `09_compare_kmer_runs.sh` | overlap of k-mer candidate peaks between runs | `09_kmer_overlap/overlap_summary.tsv` |
| 10 | `10_nucmer_dotplots.sh` | nucmer + delta-filter + show-coords + mummerplot | `10_dotplots/<pair>/` |
| 11 | `11_figures.sh` | Manhattan plots, R dotplots, chromosome correspondence | `11_figures/` |

## Standalone tools (`bin/`)

Every analysis step is also a standalone script with `--help`:

| Tool | Purpose |
|---|---|
| `classify_sex_specific_snps.py` | male-/female-specific sites from a genotype table (`--alt-dosage` for 0/0/0/1-only or 0/0/1/1-only runs) |
| `plot_sex_specific_windows.py` | two-panel genome-wide window plot + top windows |
| `merge_top_candidate_windows.py` | merge top windows into ranked candidate intervals |
| `strict_overlap_with_candidates.py` | strict sites inside relaxed candidate intervals |
| `inspect_candidate_region.py`, `run_candidate_region_inspection.py` | per-site sex bias in one / the top-N intervals |
| `plot_candidate_chromosome.py` | zoomed chromosome plot (Mb), shaded core |
| `pca_from_gt_presence.py` | PCA of a region; spot mis-sexed samples |
| `permute_sex_labels.py` | sex-label permutation p-value for a region |
| `window_enrichment_permutation.py` | is one chromosome enriched vs random window sets? |
| `plot_manhattan_windows.py` | several runs/references on one Manhattan-style plot |

`findzx/make_units.py` writes findZX sample files (XY / ZW / `--shuffle`),
and `findzx/shorten_chrom_names.sh` tidies findZX tables before replotting.

## Requirements

fastp, FastQC, MultiQC, BWA, samtools, GATK4 (read groups), freebayes,
bcftools/htslib, bedtools, KMC3, seqtk, MUMmer4 (+ gnuplot for mummerplot),
Python 3 (pandas, numpy, matplotlib, scikit-learn), R (ggplot2, dplyr), and
findZX's own environment for stage 4. All except findZX are in
`environment.yml`.

## Layout

```
run_pipeline.sh   runs the stages in order
scripts/          stage scripts 01-11 + common.sh
bin/              Python analysis and plotting tools
R/                dotplots and chromosome correspondence
findzx/           findZX helpers (units, config template, replotting)
external/findZX   findZX (git submodule)
config/           example config, sample sheet, assembly pairs
docs/WORKFLOW.md  step-by-step account of the analysis and its decisions
```

## Citing

If you use the findZX stage, cite the findZX paper as given in
`external/findZX/README.md`. Cite the other tools (BWA, freebayes, bcftools,
bedtools, KMC, MUMmer) as usual.
