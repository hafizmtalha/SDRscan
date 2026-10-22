# SDRscan validation on simulated data

SDRscan was run end to end on simulated datasets with a known sex-determining
region (SDR), and every result was checked automatically against what was
planted. The tests are in `tests/` and can be repeated on any installation
(see "Test your installation" in the README).

## The simulated genome

`tests/simulate_dataset.py` builds a 3.3 Mb genome and paired-end reads
(150 bp, with varying base qualities and sequencing errors):

| Sequence | Size | Content |
|---|---|---|
| chr_01, chr_02 | 1.2 Mb, 1.0 Mb | autosomes |
| chr_03 | 1.0 Mb | X/Y pair with an SDR at 300–700 kb |
| scaffold_1 | 60 kb | unplaced scaffold (ignored by `CHROM_REGEX`) |

Planted on the Y (XY system; males carry one Y among P haplotypes, females
none):

| Feature | Position | Expected signal |
|---|---|---|
| Y-specific SNPs, every 150 bp | chr_03:300–600 kb | male-specific SNPs, single ALT dose, higher male heterozygosity |
| Y-only 5 kb insertion | chr_03:500 kb | male-specific k-mers |
| Y deletion | chr_03:600–700 kb | lower male depth: log2 = −1 (diploid), log2(3/4) = −0.415 (tetraploid) |
| Y inversion | chr_03:350–450 kb | reverse alignment in the female-vs-male dotplot |

Autosomes carry ordinary population variation shared by both sexes, so chance
"sex-specific" variants occur as in real data.

A female reference (X only) and a male reference (Y in place of X) are
written, and the female reference is used for read mapping.

## Results

Run on Linux (Ubuntu 24.04), 4 threads, 8 GB RAM, with the settings written
by `tests/run_test.sh`. Every stage except findZX (stage 4)
was run.

| Test | Design | Checks | Run time |
|---|---|---|---|
| `diploid_2v2` | ploidy 2, 2 M + 2 F, 15× | **30/30 passed** | 6 min |
| `tetraploid_6v6` | ploidy 4, 6 M + 6 F + 2 unknown sex, 20× | **33/33 passed** | ~31 min |
| `pilot_1v1` | ploidy 2, 1 M + 1 F, 15×, `ALLOW_SINGLE_PER_SEX=1` | **14/14 passed** | 4 min |

### Diploid, 2 males + 2 females

| Stage | Result |
|---|---|
| 2–3 | joint VCF, 8,724 biallelic SNPs, 0% missing genotypes |
| 5 | 2,207 male-specific SNPs, 91% inside the SDR; top candidate interval chr_03:270–610 kb |
| 5 | 230 female-specific vs 208 male-specific SNPs outside the SDR: the chance background expected with 2 per sex |
| 6 | permutation p = 0.167, which is the floor for 2 + 2 (6 labellings): a 2 + 2 design cannot reach significance, and SDRscan says so |
| 6 | PCA separates the sexes on PC1 |
| 7 | skipped (needs ≥ 3 males), as designed |
| 8 | male-specific k-mers peak in the SDR on both references |
| 11 | chromosomes correspond 1:1; the inversion shows as 1 reverse alignment |
| 12 | lower male depth at chr_03:600–700 kb (log2 −1.00), higher male heterozygosity at 300–600 kb; nothing strong elsewhere |
| 13 | all 4 recorded sexes confirmed |

### Tetraploid, 6 males + 6 females + 2 of unknown sex

| Stage | Result |
|---|---|
| 2–3 | 14 samples called jointly at ploidy 4, 9,131 SNPs, all genotypes with 4 alleles, 0% missing |
| 5 | relaxed: 2,021 male-specific SNPs, 99% in the SDR; 19 female-specific (background) |
| 5 | strict: 2,000 male-specific (1,998 in the SDR), 1 female-specific |
| 5 | dosage classes: the 1-copy run (`0/0/0/1`) finds 2,026 sites; the 2- and 3-copy runs find none in the SDR, as expected for one Y among four copies |
| 6 | exact permutation over all 924 labellings, **p = 0.0076** (floor 0.0011) |
| 6 | presence and dosage PCAs both separate the sexes |
| 7 | the SDR is found in all 6 leave-one-male-out runs |
| 8 | male-specific k-mers peak in the SDR on both references |
| 11 | chromosomes 1:1; inversion seen |
| 12 | lower male depth at 600–700 kb, **log2 −0.41** (expected −0.415 at ploidy 4); higher male heterozygosity at 300–600 kb; no other candidates |
| 13 | both unknown samples sexed correctly with high confidence (U01 F, U02 M); all 12 recorded sexes confirmed |

### Pilot, 1 male + 1 female

| Stage | Result |
|---|---|
| 1–3 | joint VCF, 8,042 SNPs |
| 5, 6, 7, 8, 9 | skipped, as designed: without replication every private variant looks sex-specific |
| 10 | assemblies aligned |
| 12 | lower male depth at chr_03:600–700 kb (log2 −1.00) and higher male heterozygosity across 300–600 kb, both as multi-window (`strong`) runs; 47 single- or double-window `weak` hits elsewhere, as expected when comparing two individuals |
| 13 | both recorded sexes confirmed from depth alone (medium confidence) |

## Problems the tests found and fixed

Running the full pipeline on data with a known answer found these, all fixed
before release:

1. **Per-sample variant calling hid sex-specific SNPs.** Calling each sample
   separately and merging the VCFs left samples matching the reference as
   missing (`./.`) instead of 0/0. At Y-specific SNPs every female looked
   missing, the sites failed the "called in enough females" filter, and the
   2 + 2 test found **0** male-specific SNPs. Stage 2 now calls all samples
   jointly in genome chunks (see `docs/WORKFLOW.md`, "Correction: joint
   calling").
2. **Maximum thresholds were rounded up.** `--max-alt-opposite 0.25` with 2
   females allowed 1 female to carry the allele. Fractions now round towards
   the stricter side (minimums up, maximums down), and the default relaxed
   run uses a fraction (0.2) instead of a fixed count of 2, which with 2
   samples per sex would have allowed every opposite-sex sample.
3. **Sexing weighted a few chance markers like thousands of real ones.** The
   SNP and k-mer scores averaged the share of male-specific and
   female-specific markers carried. With 2,000 male markers and 1 chance
   "female" marker, a male carrying that one allele scored 0. Every marker
   now counts once; the affected sample went from a low-confidence to a
   high-confidence call.
4. **Integer overflow in the dotplot script.** Multiplying two alignment
   lengths overflowed R integers, so forward/reverse counts were NA. Fixed;
   the planted inversion is now reported.
5. **FastQC 0.13 removed the `-o` option.** Stage 1 uses `--outdir`.
6. **k-mer "present in all but one" with 2 samples meant 1 sample.** The
   new default `auto` uses `all` below 4 samples of a sex.
7. **The depth threshold and sexing expectation assumed diploidy.** Both now
   follow the declared ploidy ((P−1)/P depth in the heterogametic sex).
8. **Depth/heterozygosity candidates mixed real runs with single noisy
   windows.** Candidates now carry `support` = `strong` (≥ 3 adjacent
   windows) or `weak`.

## Not covered by these tests

- **findZX (stage 4)** is not part of these tests; it is a separately
  published workflow with its own requirements (README, "findZX (stage 4)
  requirements"). Run `RUN_FINDZX=1 tests/run_test.sh with_findzx 2 3 3` to
  include it.
- **ZW systems** were not simulated. The code paths are symmetric (female
  reference for XY, male reference for ZW, `FOCAL_SEX=F`), but only XY was
  tested.
- **Mixed-ploidy populations** (a `ploidy` column in the sample sheet,
  freebayes `--cnv-map`) were not simulated.
- **Real-data scale:** the simulated genome is 3.3 Mb. On a real plant
  genome of about 500 Mb, expect hours per sample for alignment, and
  variant-calling time that grows with genome size and sample number.
