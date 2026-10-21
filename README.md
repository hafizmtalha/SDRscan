# sex-region-pipeline

Scripts for finding the sex-determining / sex-linked region of a genome from
whole-genome resequencing of males and females. It combines three independent
lines of evidence and a set of robustness checks:

| Evidence | Idea | Reference needed? |
|---|---|---|
| **A. Depth & heterozygosity** ([findZX](https://github.com/hsigeman/findZX)) | Sex-linked regions differ between sexes in read depth and heterozygosity | yes |
| **B. Sex-specific SNP density** | Count SNPs whose ALT allele is (almost) only found in one sex, in sliding windows | yes |
| **C. Sex-specific k-mers** | k-mers present in most of one sex and absent from every individual of the other, mapped back to the genome | only for placing them |
| **D. Assembly comparison** | Whole-genome dotplots between male/female (or phased) assemblies to look for structural differences and check that chromosome numbering agrees | — |

Robustness checks: strict vs relaxed thresholds, candidate-region genotype
inspection, leave-one-out, PCA, sex-label permutation, per-chromosome
enrichment permutation, randomised-sex findZX control.

Genotypes are handled as ALT-present / ALT-absent per sample, so everything
works for diploids and polyploids (it was developed on tetraploid calls).

---

## Requirements

```bash
conda env create -f environment.yml
conda activate sex-region
```

fastp, FastQC, MultiQC, BWA, samtools, GATK4 (read groups only), freebayes,
bcftools/htslib, bedtools, KMC3, seqtk, MUMmer4, Python 3 (pandas, numpy,
matplotlib, scikit-learn), R (ggplot2, dplyr). findZX is installed separately
from its own repository.

## Setup

```bash
cp config/pipeline.env.example config/pipeline.env     # edit paths, threads, ploidy
cp config/samples.example.tsv  config/samples.tsv      # sample <tab> sex (M/F) <tab> fq1 <tab> fq2
```

All stage scripts read `config/pipeline.env` (or `CONFIG=/other/file`) and
write under `OUTDIR` (default `results/`). Run long stages with `nohup` or in
`tmux`; on shared servers keep `THREADS` and `MAX_PARALLEL_SAMPLES` modest.

---

## Workflow

```
reads ──01──> BAM ──02──> per-sample VCF ──03──> merged, filtered SNPs + genotype table
                                                     │
                     ┌───────────────────────────────┼─────────────────────────┐
                     04 sex-specific SNP windows     05 leave-one-out         bin/ PCA, permutations
reads ─────────────> 06 sex-specific k-mers ──> 07 overlap between runs
reads ─────────────> findZX (XY / ZW / shuffled)
assemblies ────────> 08 nucmer ──> R/ dotplots + chromosome correspondence
```

### 1. QC, trimming, alignment — `scripts/01_qc_align.sh`

fastp → FastQC (raw and trimmed) → `bwa mem` → `samtools sort/index` →
GATK `AddOrReplaceReadGroups` (sample name = sample ID) → MultiQC.
Set `MARK_DUPLICATES=1` to also run `MarkDuplicatesSpark`.

### 2. Variant calling — `scripts/02_call_variants.sh`

`freebayes -p PLOIDY -C 4 --min-mapping-quality 20 --min-base-quality 20
--genotype-qualities --use-mapping-quality`, one sample per job,
`MAX_PARALLEL_SAMPLES` at a time. Use `PLOIDY=4` if the population may be
tetraploid. Empty VCFs are reported as errors (a truncated BAM gave 0-byte
output for one sample).

### 3. Merge and filter — `scripts/03_merge_vcfs.sh`

`bcftools norm -m -both` → rename the sample to its ID → `bcftools merge` →
`QUAL >= 30 && INFO/DP >= 5` → SNPs only → biallelic only. Writes:

- `merged.filtered.snps.vcf.gz`
- `merged.genotypes.tsv` (`bcftools query -f '%CHROM\t%POS[\t%GT]\n'`)
- `vcf_samples.txt` (column order), `males.txt`, `females.txt`

### 4. Sex-specific SNP density — `scripts/04_sex_specific_snp_windows.sh RUN [options]`

1. **Classify** (`bin/classify_sex_specific_snps.py`). A site is
   *male-specific* if ≥ `--min-alt-males` males carry ALT and ≤
   `--max-alt-opposite` females do (both sexes having enough called
   genotypes), and vice versa for *female-specific*.
2. Sites → BED → `bedtools makewindows` (`WINDOW`/`STEP`, main chromosomes
   only) → `bedtools intersect -c`.
3. Top 50 windows per sex → merged candidate intervals
   (`bedtools merge` and `bin/merge_top_candidate_windows.py`, which reports
   `N_WINDOWS`, `MAX_COUNT`, `MEAN_COUNT`) → counts per chromosome →
   genome-wide plot.

Run it several times with different thresholds (one `RUN` name each). How
strict to be depends on sample size — aim for a few hundred sites, not 4 and
not 50,000:

```bash
# relaxed: discovery
scripts/04_sex_specific_snp_windows.sh relaxed \
  --min-called-males 6 --min-called-females 5 --min-alt-males 6 --min-alt-females 5 --max-alt-opposite 2
# strict: validation (no ALT allowed in the other sex)
scripts/04_sex_specific_snp_windows.sh strict \
  --min-called-males 7 --min-called-females 5 --min-alt-males 7 --min-alt-females 5 --max-alt-opposite 0
# genotype-class runs (polyploids): count only single-dose or only double-dose ALT genotypes
scripts/04_sex_specific_snp_windows.sh dose1 ... --alt-dosage 1     # 0/0/0/1
scripts/04_sex_specific_snp_windows.sh dose2 ... --alt-dosage 2     # 0/0/1/1
```

To test an XY system, call against a reference from the homogametic sex
(female) and look for **male**-specific sites; for ZW use a male reference and
look for **female**-specific sites. Running both is a useful contrast.

#### Following up candidates

| Script | Question it answers |
|---|---|
| `bin/strict_overlap_with_candidates.py` | Do strict sites fall inside the relaxed-threshold candidate intervals? |
| `bin/inspect_candidate_region.py` / `bin/run_candidate_region_inspection.py` | Is a peak made of many repeated sex-biased SNPs or a few noisy sites? |
| `bin/plot_candidate_chromosome.py` | What does the signal look like along one chromosome (Mb axis, shaded core)? |
| `bin/pca_from_gt_presence.py` | Do samples split by sex in the candidate region? Does any sample cluster with the other sex (mis-sexing)? |
| `bin/permute_sex_labels.py` | How often do random sex labellings give as many strong sex-biased sites in the region? |
| `bin/window_enrichment_permutation.py` | Is one chromosome enriched for sex-specific SNPs/k-mers compared with random window sets? |
| `scripts/05_leave_one_out.sh` | Does the candidate survive dropping each sample of the focal sex? |
| `bin/plot_manhattan_windows.py` | Several runs/references on one Manhattan-style figure |

Example: subset a region and run the PCA

```bash
cd results/03_merged
bcftools view -r chr_07:8000000-9000000 merged.filtered.snps.vcf.gz -Oz -o region.vcf.gz
bcftools query -f '%CHROM\t%POS[\t%GT]\n' region.vcf.gz > region.genotypes.tsv
python ../../bin/pca_from_gt_presence.py --genotypes region.genotypes.tsv --outprefix region_PCA \
  --min-called 12 --min-sites 100 --max-sites 20000
```

### 5. Sex-specific k-mers — `scripts/06_kmer_sex_specific.sh MSK|FSK MIN_PRESENT TAG=ref.fa ...`

1. KMC (`k=31`, `-ci2`) per sample.
2. Focal-sex set: k-mers present in ≥ `MIN_PRESENT` focal individuals.
3. Subtract the union of all opposite-sex k-mers → **MSK** (male-specific) or
   **FSK** (female-specific) k-mers.
4. Subsample 200k, `bwa mem` to each reference, keep MAPQ ≥ 20.
5. Count per 25 kb / 5 kb window and per chromosome; merge the top 200
   windows into candidate peaks.

Map to every available assembly/haplotype: a sex-limited chromosome (Y/W) may
only be assembled in one of them.

```bash
scripts/06_kmer_sex_specific.sh MSK 7 F1=female_h1.fa F2=female_h2.fa M1=male_h1.fa M2=male_h2.fa
scripts/06_kmer_sex_specific.sh FSK 6 F1=female_h1.fa F2=female_h2.fa M1=male_h1.fa M2=male_h2.fa
scripts/07_compare_kmer_runs.sh results/07_kmer_overlap \
  MSK_F1=results/06_kmer/MSK/windows/kmer_MSK_candidates.F1.bed \
  FSK_F1=results/06_kmer/FSK/windows/kmer_FSK_candidates.F1.bed
```

### 6. findZX — `findzx/`

```bash
python findzx/make_units.py --samples config/samples.tsv --system XY --out <findZX>/config/units_XY.tsv
python findzx/make_units.py --samples config/samples.tsv --system ZW --out <findZX>/config/units_ZW.tsv
python findzx/make_units.py --samples config/samples.tsv --system XY --shuffle --seed 7 \
  --out <findZX>/config/units_shuffled.tsv       # negative control
cd <findZX> && /path/to/sex-region-pipeline/findzx/run_findzx.sh config/my_config.yml 16
```

`findzx/config.example.yml` has the settings we used. Pair XY runs with a
female reference and ZW runs with a male reference. After a run, check that
the plots contain every sample. `findzx/shorten_chrom_names.sh` shortens
chromosome names in the output tables before replotting.

### 7. Assembly comparison — `scripts/08_nucmer_dotplots.sh pairs.tsv`

`nucmer --mum` → `delta-filter -i 95 -l 1000 -q` → `show-coords` →
`mummerplot`, for every pair in `config/pairs.example.tsv` (set
`CHROM_REGEX` to drop unplaced scaffolds). Then:

```bash
Rscript R/dotplot_genome.R results/08_dotplots/A_vs_B/A_vs_B.filtered95.coords A B A_vs_B_genome
Rscript R/dotplot_chromosome.R results/08_dotplots/A_vs_B/A_vs_B.filtered95.coords \
  A_chr_08 B_chr_08 A B A_vs_B_chr8 10000 5
```

`dotplot_genome.R` also writes a **chromosome correspondence table** (which
query chromosome receives most of each reference chromosome's aligned
sequence). Check it before comparing results across assemblies, because
chromosome numbers are not guaranteed to match between assemblies.

---

## Interpreting results

Treat a region as a strong candidate only when several independent lines of
evidence agree. Sex-specific SNP density, sex-specific k-mers and findZX
depth/heterozygosity are affected differently by sample size, repeats and
mapping artefacts, so agreement between them carries more weight than any
single peak. With small samples (around 10 individuals), single windows with
2–3 sex-specific sites can be chance; use the strict overlap, leave-one-out
and permutation outputs to judge them.

## Repository layout

```
bin/        Python tools (classification, windows, follow-up, statistics, plots)
scripts/    Stage scripts 01–08 and shared common.sh
findzx/     Helpers for running findZX (units, config, replotting)
R/          Dotplot and chromosome correspondence plots
config/     Example config, sample sheet and assembly pair list
```
