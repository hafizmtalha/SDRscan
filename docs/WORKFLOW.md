# SDRscan workflow: how the analysis was built, step by step

This document follows the analysis in the order it was actually done. For
each step it gives the question, the script that answers it, the key
settings, and what the result led us to do next. The aim is that someone can
repeat the same reasoning on a new dataset, not only rerun the commands.

The study design: whole-genome resequencing of wild males and females (about
10 individuals at first, 16 later), several reference assemblies (two phased
haplotypes of a female genome, two of a male genome, one more of unknown
sex), and a species whose populations can be diploid or tetraploid. The
question was which chromosome region carries the sex-determining region
(SDR), and whether the system is XY or ZW.

```
 reads
   │  Step 1  QC, trimming, alignment ................ 01_qc_align.sh
   ▼
 BAM ── Step 2  variant calling (diploid → tetraploid) 02_call_variants.sh
   │
   ├──── Step 3  findZX: depth + heterozygosity ......... 04_findzx.sh  ──────────┐
   ▼                                                                             │
 VCF ── Step 4  merge + filter ......................... 03_merge_vcfs.sh         │
   │                                                                             │
   │    Step 5  sex-specific SNP density ............... 05_sex_specific_snp_windows.sh
   │    Step 6  is the peak real? ...................... 06_candidate_followup.sh, 07_leave_one_out.sh
   │                                                                             │
 reads ─ Step 7  sex-specific k-mers .................... 08_kmer_sex_specific.sh, 09_compare_kmer_runs.sh
                                                                                 │
 assemblies ─ Step 8  structure + chromosome naming ..... 10_nucmer_dotplots.sh, R/
                                                                                 │
 Step 9  more samples, both references, both systems ───────────────────────────┤
 Step 10 final figures ................................. 11_figures.sh ◄─────────┘
```

---

## Step 1 — Read QC and alignment

**Question:** are the reads usable, and where do they map?

**Script:** `scripts/01_qc_align.sh` (stage 1)

fastp trimming → FastQC before and after → `bwa mem` to the reference →
`samtools sort/index` → GATK `AddOrReplaceReadGroups`, with the sample ID as
read group and sample name → MultiQC report.

**Settings and choices**
- The full reference assembly was used for alignment, including unplaced
  scaffolds, so reads from them are not forced onto the chromosomes. The
  scaffolds are filtered out later at the window/plot stage (`CHROM_REGEX`).
- Read groups matter: without them every VCF sample is called `unknown`, and
  `bcftools merge` then fails or merges samples into one.
- We did not mark duplicates or recalibrate base qualities in the runs used
  for results (`MARK_DUPLICATES=0`). State this in the methods.

**Led to:** BAMs for every sample → variant calling.

---

## Step 2 — Variant calling, first diploid, then tetraploid

**Question:** what are each sample's genotypes?

**Script:** `scripts/02_call_variants.sh` (stage 2)

```
freebayes -f REF -p PLOIDY -C 4 --min-mapping-quality 20 --min-base-quality 20 \
          --genotype-qualities --use-mapping-quality sample.RG.bam
```

**What happened**
1. A single joint run (`freebayes-parallel` over all BAMs, 60 threads) took
   down the shared server. We switched to **one sample per job, 2–3 jobs at a
   time** (`MAX_PARALLEL_SAMPLES`), which is what the script does.
2. The first calls were diploid (`-p 2`). The population could be
   tetraploid, so **every sample was re-called with `-p 4`**, and the
   tetraploid calls were used from then on.
3. Because of (2), ploidy is now a required setting (`PLOIDY`, or a
   per-sample `ploidy` column) with no default. Stage 2 records the ploidy
   each sample was called at, stage 3 refuses to merge if the config has
   changed since, and every later tool reads genotypes at that ploidy. A
   genotype with the wrong number of alleles is reported instead of being
   silently miscounted.
4. One sample repeatedly produced an empty VCF. The cause was its BAM;
   re-aligning fixed it. The script now stops on a 0-byte VCF instead of
   letting it through silently.

**Led to:** per-sample tetraploid VCFs → merging (Step 4). findZX (Step 3)
ran in parallel because it works from the raw reads.

---

## Step 3 — findZX: depth and heterozygosity

**Question:** does any region differ between the sexes in read depth (one
sex has one copy, so half the depth) or heterozygosity?

**Script:** `scripts/04_findzx.sh` (stage 4) runs the
[findZX](https://github.com/hsigeman/findZX) Snakemake pipeline from the
`external/findZX` submodule. `findzx/make_units.py` builds findZX's sample
file, and `findzx/config.example.yml` holds the settings.

**Runs**

| Run | Reference | System | Why |
|---|---|---|---|
| XY | female assembly | males heterogametic | main hypothesis |
| ZW | male assembly | females heterogametic | the alternative system |
| SHUFFLED | female assembly | XY with random sex labels | negative control: the real pattern should disappear |

**Lessons**
- It took several attempts to get findZX running on a current system: old
  pinned tools clashing with newer libraries, a dependency that didn't
  install, a config syntax error. The fixes that made it work:
  - a dedicated `--conda-prefix`
  - `--conda-frontend mamba`
  - `--rerun-incomplete`
  - re-running after each fix (Snakemake resumes)

  The test dataset that comes with findZX is a quick way to check the
  installation before using real data.
- Pair the system with the right reference: test XY on a female reference
  and ZW on a male reference.
- **Check that the output plots show every sample.** Two runs silently used
  an older 10-sample units file after more samples had been added.
- Plots are easier to read with only assembled chromosomes (`chr_file`) and
  short chromosome names (`findzx/shorten_chrom_names.sh`, then rerun only
  the plotting rule).
- A findZX run on 16 samples took about two weeks on 16 cores, so start it
  early and run it in the background (`FINDZX_BACKGROUND=1`).

**Led to:** findZX gives a broad picture. To pin down a region at SNP
resolution we built the sex-specific SNP analysis (Steps 4–5).

---

## Step 4 — Merge and filter the VCFs

**Question:** one clean multi-sample SNP table to work from.

**Script:** `scripts/03_merge_vcfs.sh` (stage 3)

1. `bcftools norm -m -both` per sample and rename the sample to its ID
   (`bcftools reheader`)
2. `bcftools merge` all samples
3. `QUAL >= 30 && INFO/DP >= 5`
4. SNPs only (`-v snps`), biallelic only (`-m2 -M2`)
5. Export genotypes:
   ```
   bcftools query -f '%CHROM\t%POS[\t%GT]\n' merged.filtered.snps.vcf.gz > merged.genotypes.tsv
   bcftools query -l merged.filtered.snps.vcf.gz > vcf_samples.txt
   ```

The genotype table, `vcf_samples.txt` (column order) and `males.txt` /
`females.txt` are the inputs to everything SNP-based.

**Why not FST?** Windowed Weir & Cockerham FST (vcftools) assumes diploid
genotypes, so it was dropped once we moved to tetraploid calls. The
presence/absence approach below works at any ploidy.

---

## Step 5 — Sex-specific SNP density

**Question:** where do SNPs cluster whose ALT allele is carried by (almost)
all of one sex and (almost) none of the other? In an XY system, Y-linked
variants are male-specific. In ZW, W-linked variants are female-specific.

**Script:** `scripts/05_sex_specific_snp_windows.sh RUN [classifier options]`
(stage 5, one run per `SNP_RUNS` entry)

1. **Classify** each site (`bin/classify_sex_specific_snps.py`). Each
   genotype is reduced to ALT-present / ALT-absent / missing, so ploidy does
   not matter. A site is *male-specific* when:
   - at least `min_called` males and females have a genotype,
   - at least `min_alt_males` males carry ALT,
   - at most `max_alt_opposite` females carry ALT.

   *Female-specific* is the reverse.
2. Sites → BED → `bedtools makewindows -w 100000 -s 50000` →
   `bedtools intersect -c`.
3. Top 50 windows per sex → merge overlapping ones into candidate intervals
   (`bin/merge_top_candidate_windows.py` reports the number of windows, max
   and mean count; `bedtools merge` gives the plain intervals) → counts per
   chromosome.
4. Genome-wide plot (`bin/plot_sex_specific_windows.py`).

**Choosing thresholds.** This was the main tuning step:

| Setting | Sites | Verdict |
|---|---|---|
| all of one sex ALT, none of the other | 4 of 23 million | too strict, nothing to plot |
| most of one sex, ≤ 2 of the other ("relaxed") | ~1,500 | good for **finding** candidates |
| most of one sex, 0 of the other ("strict") | ~80 | good for **validating** them |

Aim for a few hundred sites. Scale the thresholds with sample size: with
9 males and 7 females we used about 7/5 for the focal/opposite sex. The
classifier now accepts fractions (`--min-alt-males 0.8`), so a threshold set
carries over to a study of a different size. `SNP_RUNS` in the config holds
any number of named threshold sets, and each run logs the counts the
fractions resolved to.

**Genotype-class runs.** With tetraploid calls, "ALT present" can mean one
ALT copy (0/0/0/1) or two (0/0/1/1), which behave differently.
`--alt-dosage 1` and `--alt-dosage 2` count only one class each, which shows
whether the signal comes from single-dose variants, as expected for a
sex-limited haplotype. For any ploidy P, `SNP_DOSAGE_CLASSES=auto` runs every
class from 1 to P−1.

**Led to:** a ranked list of candidate intervals on several chromosomes, each
supported by only a handful of windows → is any of them real? (Step 6)

---

## Step 6 — Is the peak real?

**Question:** a window with 3–5 sex-specific SNPs can be chance. Which
candidates survive scrutiny?

**Scripts:** `scripts/06_candidate_followup.sh` (stage 6, regions from
`CANDIDATES`) and `scripts/07_leave_one_out.sh` (stage 7). Each check is
also a standalone tool in `bin/`.

| Check | Tool | What a real SDR looks like |
|---|---|---|
| Strict sites inside relaxed intervals | `strict_overlap_with_candidates.py` | strict sex-specific sites fall in the interval |
| Zoomed density | `plot_candidate_chromosome.py` | one contiguous block, not scattered spikes |
| Genotype inspection | `inspect_candidate_region.py` | many repeated sex-biased SNPs, not 1–2 noisy sites |
| PCA of the region | `pca_from_gt_presence.py` | samples split by sex; a sample with the other sex suggests mis-sexing |
| Sex-label permutation | `permute_sex_labels.py` | real labels give more strong sites than random labels (small p) |
| Chromosome enrichment | `window_enrichment_permutation.py` | candidate chromosome has more sex-specific windows than random window sets |
| Leave-one-out | `07_leave_one_out.sh` | candidate stays when each focal-sex sample is dropped |

**What we found:** one interval of a few hundred kb kept coming back. It
appeared in the relaxed runs, held a strict site, and had a block of strong
male-biased SNPs. With only about 10 samples, the conclusion was that it
*could still be noise*. That called for evidence that doesn't depend on
SNP calling, and for more samples.

---

## Step 7 — Sex-specific k-mers

**Question:** is there sequence present in one sex and absent in the other?
This needs no SNP calling, and less of a reference (only to place the
k-mers).

**Scripts:** `scripts/08_kmer_sex_specific.sh MSK|FSK MIN_PRESENT TAG=ref.fa ...`
(stage 8) and `scripts/09_compare_kmer_runs.sh` (stage 9).

1. KMC per sample: k = 31, minimum count 2 (drops most sequencing errors).
2. Focal set: k-mers present in ≥ `MIN_PRESENT` individuals of the focal sex
   (all but one, to tolerate one low-coverage or mis-sexed sample).
3. Subtract the union of every opposite-sex k-mer → **MSK** (male-specific)
   or **FSK** (female-specific).
4. Subsample 200,000, map with `bwa mem` to each reference, keep MAPQ ≥ 20.
5. Count per 25 kb / 5 kb window and per chromosome; merge the top windows
   into candidate peaks.

**Runs:** MSK and FSK, each mapped to all four haplotypes (8 runs). The
haplotypes matter: a Y or W haplotype may only be assembled in one of them,
so male-specific k-mers mapped to a female assembly can land in odd places.
Stage 9 counts overlapping candidate peaks between runs. Chromosome names are
normalised first so the two haplotypes of the same assembly can be compared.

**What we found:** the first MSK run (one reference) was scattered and
inconclusive. FSK runs and the full set across references were more
informative and pointed to a different chromosome than the SNP analysis. The
per-chromosome counts and enrichment permutation (Step 6 tool) were used to
test this.

---

## Step 8 — Assembly structure and chromosome naming

**Questions:** do the male and female assemblies differ structurally on the
candidate chromosome (inversions, rearrangements, missing blocks)? And are
"chromosome 8" in one assembly and "chromosome 8" in another really the same
chromosome?

**Scripts:** `scripts/10_nucmer_dotplots.sh pairs.tsv` (stage 10), then
`R/dotplot_genome.R` and `R/dotplot_chromosome.R` (stage 11).

1. Keep assembled chromosomes only (`CHROM_REGEX`).
2. `nucmer --mum --mincluster 100 --minmatch 20`
3. `delta-filter -i 95 -l 1000 -q`. 95% identity was a good balance; 98% was
   also tried.
4. `show-coords -T -r -c -l` → R plots:
   - genome-wide dotplot with Chr1..ChrN labels and a **chromosome
     correspondence table** (which query chromosome gets most of each
     reference chromosome's aligned bases)
   - single-chromosome dotplots with Mb axes for all pairs of the candidate
     chromosome

**Lessons**
- **Always build the correspondence table before comparing results across
  assemblies.** Here, the published SDR chromosome and our SNP candidate had
  different numbers, and the first check was whether the numbering differed.
  The table showed the numbering matched (and caught one wrongly downloaded
  assembly).
- minimap2 (asm5 + PAF dotplots) gave fragmented, low-coverage alignments
  for these assemblies; nucmer was clearer.
- For single-chromosome figures, points (not segments) and a minimum
  alignment length removed most visual noise.
- Synteny between assemblies at gene level (riparian plots) was done with
  GENESPACE, from the assemblies' protein and annotation files. That is
  outside this repository.

---

## Step 9 — More samples, both references, both systems

With the first pass done, the analysis was repeated on **16 samples (9 male,
7 female)**:

1. New samples went through Steps 1–2. Variant calling was done twice, on
   a female reference and on a male reference, giving two
   `OUTDIR`s/configs.
2. Step 4 merge, then Step 5 with thresholds rescaled to 16 samples.
3. Female reference → look for **male-specific** sites (XY test). Male
   reference → **female-specific** sites (ZW test).
4. Genotype-class runs (0/0/0/1 only, 0/0/1/1 only) on both, with 10 kb
   windows for finer resolution.
5. findZX: XY on the female reference, ZW on the male reference, and the
   shuffled-sex control (Step 3).

```bash
CONFIG=config/female_ref.env ./run_pipeline.sh 1 2 3 5
CONFIG=config/male_ref.env   ./run_pipeline.sh 1 2 3 5
```

---

## Step 10 — Final figures

**Script:** `scripts/11_figures.sh` (stage 11) and
`bin/plot_manhattan_windows.py` directly for combined panels.

- Manhattan-style SNP density: male- and female-specific windows, both
  references on one panel (colour = sex, shape = reference)
- Manhattan-style k-mer density: MSK on all four haplotypes on one plot, and
  FSK the same
- findZX genome-wide and per-chromosome plots (from findZX itself)
- Genome-wide and candidate-chromosome dotplots

Example for a combined SNP figure across two references:

```bash
python bin/plot_manhattan_windows.py \
  --track "Ref F male-specific SNPs=results_F/05_snp_windows/dose1/male_specific_windows_10000_5000.bed:green:o" \
  --track "Ref F female-specific SNPs=results_F/05_snp_windows/dose1/female_specific_windows_10000_5000.bed:purple:o" \
  --track "Ref M male-specific SNPs=results_M/05_snp_windows/dose1/male_specific_windows_10000_5000.bed:blue:^" \
  --track "Ref M female-specific SNPs=results_M/05_snp_windows/dose1/female_specific_windows_10000_5000.bed:orange:^" \
  --chrom-regex 'chr_?[0-9]+$' --title "Genome-wide sex-specific SNP density (0/0/0/1)" \
  --outprefix snp_dose1_both_refs
```

---

## Things we would do differently

- **Log every command as it is run.** Most of the time spent reconstructing
  this workflow went into recovering command order. The stage scripts write
  their options and logs into each output folder for this reason.
- **Fix the thresholds per sample size before looking at results**, and
  report every threshold set used.
- **Run the negative controls early:** shuffled-sex findZX and the label
  permutation. They are cheap compared with chasing a chance peak.
- **Mark duplicates.** It wasn't done in the original runs, and it reduces
  false heterozygous calls.
- **Size jobs for the server.** Several crashes and restarts came from
  launching too many threads or too many findZX runs at once.

---

## Added after the original analysis: pilot runs and sex assignment

The original study always had several individuals of each sex. Two modes
were added later so SDRscan is useful at both ends of the sample-size range.
They were tested on simulated data, not on the original study.

**Pilot (1 male + 1 female), `ALLOW_SINGLE_PER_SEX=1`.** Only analyses that
make sense without replication run: findZX and the new depth/heterozygosity
scan (stage 12). The scan normalises each sample's depth by its own genome
median, reports log2(male/female) per window, and marks runs of windows with
a strong shift. A well-differentiated X or Z shows about half depth in the
heterogametic sex even in one individual. SNP density, k-mers and every test
are skipped, because with one individual per sex every private variant looks
sex-specific. Stage 12 also runs for larger designs, where it requires every
male to differ from every female in a window. Unlike findZX's heterozygosity,
it respects the declared ploidy.

**Sex assignment (stage 13).** Once a sex-linked region or sex-specific
markers are known, from this project or an earlier one, samples marked `U`
can be sexed from depth over the region, genotypes at the marker SNPs and
presence of sex-specific k-mers. Scores are calibrated on the known-sex
samples. The same step re-checks every recorded sex and flags mismatches.
This would have been useful in the original study, where possible mis-sexing
was raised to explain weak signals.
