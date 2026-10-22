#!/usr/bin/env python3
"""
Simulate a small XY test dataset for SDRscan: reference assemblies, paired-end
reads for males, females and unknown-sex samples at any ploidy, a sample
sheet, a config file and an assembly pair list.

Genome (~3.3 Mb):
  chr_01, chr_02   autosomes with population SNPs
  chr_03           sex chromosome; sex-determining region (SDR) at 300-700 kb
  scaffold_1       small unplaced scaffold (should be ignored by windows/plots)

The Y haplotype differs from X only inside the SDR:
  - fixed Y-specific SNPs (1 per ~150 bp)        -> sex-specific SNPs, dosage 1
  - a 5 kb Y-only insertion at 500 kb            -> male-specific k-mers
  - a deletion of 600-700 kb on Y                -> lower depth in males there
  - an inversion of 350-450 kb on Y              -> visible in assembly dotplots
Outside the SDR, X and Y are interchangeable (pseudo-autosomal).

A male carries one Y and (ploidy - 1) X haplotypes; a female carries
ploidy X haplotypes. Two assemblies are written: female.fa (X version of
chr_03) and male.fa (Y version of chr_03, renamed chr_03).

Usage:
  simulate_dataset.py --outdir test_run --ploidy 4 --males 6 --females 6 --unknown 2 --coverage 20
"""

import argparse
import gzip
import os
import random

BASES = "ACGT"
COMP = str.maketrans("ACGT", "TGCA")
SDR = (300_000, 700_000)
Y_INSERT_AT, Y_INSERT_LEN = 500_000, 5_000
Y_DELETION = (600_000, 700_000)
Y_INVERSION = (350_000, 450_000)


def rand_seq(n, rng):
    return "".join(rng.choices(BASES, k=n))


def mutate(base, rng):
    return rng.choice([b for b in BASES if b != base])


def revcomp(s):
    return s.translate(COMP)[::-1]


def write_fasta(path, seqs):
    with open(path, "w") as f:
        for name, s in seqs.items():
            f.write(f">{name}\n")
            for i in range(0, len(s), 80):
                f.write(s[i:i + 80] + "\n")


def build_population(rng, sizes, snp_every=400):
    """Reference sequences plus polymorphic sites with an allele frequency each."""
    ref = {c: rand_seq(n, rng) for c, n in sizes.items()}
    sites = {}
    for c, s in ref.items():
        pos = sorted(rng.sample(range(len(s)), len(s) // snp_every))
        sites[c] = [(p, mutate(s[p], rng), rng.uniform(0.1, 0.9)) for p in pos]
    return ref, sites


def haplotype(ref, sites, rng):
    """One random haplotype drawn from the population allele frequencies."""
    h = {}
    for c, s in ref.items():
        seq = list(s)
        for p, alt, freq in sites[c]:
            if rng.random() < freq:
                seq[p] = alt
        h[c] = "".join(seq)
    return h


def to_y(x_chr3, y_snps, insert):
    """Turn an X copy of chr_03 into a Y copy: SDR SNPs, inversion, insertion, deletion."""
    s = list(x_chr3)
    for p, b in y_snps:
        s[p] = b
    s = "".join(s)
    a, b = Y_INVERSION
    s = s[:a] + revcomp(s[a:b]) + s[b:]
    d0, d1 = Y_DELETION
    s = s[:d0] + s[d1:]
    return s[:Y_INSERT_AT] + insert + s[Y_INSERT_AT:]


def simulate_reads(haps, coverage, out1, out2, rng, read_len=150, insert=350, sd=40, err=0.002):
    total = sum(len(s) for h in haps for s in h.values())
    n_pairs = int(coverage * total / (2 * read_len) / len(haps)) * len(haps)
    seqs = [s for h in haps for s in h.values()]
    weights = [len(s) for s in seqs]
    # realistic-looking Phred+33 qualities (mostly Q37, some lower)
    quals = ["".join(rng.choice("FFFFFF:,") for _ in range(read_len)) for _ in range(64)]

    def add_errors(r):
        if rng.random() > err * read_len * 2:
            return r
        r = list(r)
        for _ in range(max(1, int(rng.expovariate(1 / (err * read_len))))):
            i = rng.randrange(read_len)
            r[i] = mutate(r[i], rng) if r[i] in BASES else "A"
        return "".join(r)

    with gzip.open(out1, "wt", compresslevel=3) as f1, gzip.open(out2, "wt", compresslevel=3) as f2:
        for i in range(n_pairs):
            s = rng.choices(seqs, weights)[0]
            frag = max(read_len + 10, int(rng.gauss(insert, sd)))
            if len(s) <= frag:
                continue
            start = rng.randrange(len(s) - frag)
            fr = s[start:start + frag]
            if rng.random() < 0.5:
                fr = revcomp(fr)
            r1, r2 = add_errors(fr[:read_len]), add_errors(revcomp(fr)[:read_len])
            f1.write(f"@r{i}/1\n{r1}\n+\n{quals[i % 64]}\n")
            f2.write(f"@r{i}/2\n{r2}\n+\n{quals[(i + 7) % 64]}\n")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--ploidy", type=int, default=2)
    ap.add_argument("--males", type=int, default=2)
    ap.add_argument("--females", type=int, default=2)
    ap.add_argument("--unknown", type=int, default=0, help="extra samples written with sex U (truth in truth.tsv)")
    ap.add_argument("--coverage", type=float, default=15, help="total read depth per individual")
    ap.add_argument("--seed", type=int, default=1)
    args = ap.parse_args()

    rng = random.Random(args.seed)
    out = os.path.abspath(args.outdir)
    os.makedirs(f"{out}/reads", exist_ok=True)
    os.makedirs(f"{out}/refs", exist_ok=True)

    sizes = {"chr_01": 1_200_000, "chr_02": 1_000_000, "chr_03": 1_000_000, "scaffold_1": 60_000}
    ref, sites = build_population(rng, sizes)
    # population SNPs inside the SDR would blur Y-specific signal; keep the SDR clean on X
    sites["chr_03"] = [s for s in sites["chr_03"] if not SDR[0] <= s[0] < SDR[1]]
    y_snps = [(p, mutate(ref["chr_03"][p], rng)) for p in range(SDR[0], SDR[1], 150)]
    y_insert = rand_seq(Y_INSERT_LEN, rng)

    write_fasta(f"{out}/refs/female.fa", ref)
    male_ref = dict(ref)
    male_ref["chr_03"] = to_y(ref["chr_03"], y_snps, y_insert)
    write_fasta(f"{out}/refs/male.fa", male_ref)

    plan = ([("M", "M")] * args.males + [("F", "F")] * args.females
            + [("U", rng.choice("MF")) for _ in range(args.unknown)])
    rows, truth = [], []
    n = {"M": 0, "F": 0, "U": 0}
    for recorded, true_sex in plan:
        n[recorded] += 1
        name = f"{recorded}{n[recorded]:02d}"
        haps = [haplotype(ref, sites, rng) for _ in range(args.ploidy)]
        if true_sex == "M":
            haps[0]["chr_03"] = to_y(haps[0]["chr_03"], y_snps, y_insert)
        r1, r2 = f"{out}/reads/{name}_R1.fastq.gz", f"{out}/reads/{name}_R2.fastq.gz"
        simulate_reads(haps, args.coverage, r1, r2, rng)
        rows.append(f"{name}\t{recorded}\t{r1}\t{r2}\tNA")
        truth.append(f"{name}\t{true_sex}")
        print(f"simulated {name} (true sex {true_sex}, ploidy {args.ploidy})")

    with open(f"{out}/samples.tsv", "w") as f:
        f.write("sample\tsex\tfq1\tfq2\tploidy\n" + "\n".join(rows) + "\n")
    with open(f"{out}/truth.tsv", "w") as f:
        f.write("sample\ttrue_sex\n" + "\n".join(truth) + "\n")
    with open(f"{out}/pairs.tsv", "w") as f:
        f.write(f"{out}/refs/female.fa\t{out}/refs/male.fa\tfemale_vs_male\n")
    with open(f"{out}/expected.txt", "w") as f:
        f.write(f"system\tXY\nsdr\tchr_03:{SDR[0]}-{SDR[1]}\n"
                f"male_low_depth\tchr_03:{Y_DELETION[0]}-{Y_DELETION[1]}\n"
                f"inversion\tchr_03:{Y_INVERSION[0]}-{Y_INVERSION[1]}\nploidy\t{args.ploidy}\n")
    print(f"Dataset written to {out}")


if __name__ == "__main__":
    main()
