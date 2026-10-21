#!/usr/bin/env python3
"""
PCA of a region, coloured by recorded sex, for any ploidy.

--encoding presence (default): 1 if the genotype carries any ALT allele, else 0.
--encoding dosage: ALT copies / ploidy (0..1), which keeps the extra
information polyploid genotypes carry (0/0/0/1 vs 0/0/1/1).

Sites with fewer than --min-called called samples are dropped; remaining
missing values are imputed with the site mean. Sex labels are not used to fit
the PCA, so a sample clustering with the other sex hints at mis-sexing or at
no sex-linked structure in the region.
"""

import argparse
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA

from gtutils import GenotypeParser, load_ploidy, read_list, threshold


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--genotypes", required=True)
    ap.add_argument("--vcf-samples", default="vcf_samples.txt")
    ap.add_argument("--males", default="males.txt")
    ap.add_argument("--females", default="females.txt")
    ap.add_argument("--outprefix", default="region_PCA")
    pl = ap.add_mutually_exclusive_group(required=True)
    pl.add_argument("--ploidy", type=int)
    pl.add_argument("--sample-ploidy")
    ap.add_argument("--encoding", choices=["presence", "dosage"], default="presence")
    ap.add_argument("--min-called", default="0.8", help="count or fraction of all samples (default 0.8)")
    ap.add_argument("--min-sites", type=int, default=100)
    ap.add_argument("--max-sites", type=int, default=None, help="randomly subsample to at most this many sites")
    ap.add_argument("--seed", type=int, default=1)
    args = ap.parse_args()

    samples = read_list(args.vcf_samples)
    males = set(read_list(args.males))
    females = set(read_list(args.females))
    ploidy = load_ploidy(samples, args.ploidy, args.sample_ploidy)
    parser = GenotypeParser(ploidy)
    min_called = threshold(args.min_called, len(samples))

    rows = []
    with open(args.genotypes) as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            vals = []
            for s, gt in zip(samples, parts[2:]):
                d = parser.dosage(gt, s)
                if d is None:
                    vals.append(np.nan)
                elif args.encoding == "presence":
                    vals.append(float(d > 0))
                else:
                    vals.append(d / ploidy[s])
            v = np.array(vals)
            if np.sum(~np.isnan(v)) < min_called:
                continue
            v[np.isnan(v)] = np.nanmean(v)
            rows.append(v)

    X = np.array(rows, dtype=float)
    if args.max_sites and X.shape[0] > args.max_sites:
        idx = np.random.default_rng(args.seed).choice(X.shape[0], args.max_sites, replace=False)
        X = X[np.sort(idx)]
    if X.shape[0] < args.min_sites:
        parser.report(sys.stderr)
        raise SystemExit(f"Too few sites after filtering: {X.shape[0]}")

    pca = PCA(n_components=2)
    pcs = pca.fit_transform(X.T)

    out = pd.DataFrame({
        "sample": samples,
        "ploidy": [ploidy[s] for s in samples],
        "PC1": pcs[:, 0],
        "PC2": pcs[:, 1],
        "label": ["male" if s in males else "female" if s in females else "unknown" for s in samples],
    })
    out.to_csv(f"{args.outprefix}.scores.tsv", sep="\t", index=False)

    plt.figure(figsize=(7, 6))
    for lab in ("male", "female", "unknown"):
        sub = out[out["label"] == lab]
        if len(sub):
            plt.scatter(sub["PC1"], sub["PC2"], label=lab)
    for _, r in out.iterrows():
        plt.text(r["PC1"], r["PC2"], r["sample"], fontsize=8)
    plt.xlabel(f"PC1 ({pca.explained_variance_ratio_[0] * 100:.1f}%)")
    plt.ylabel(f"PC2 ({pca.explained_variance_ratio_[1] * 100:.1f}%)")
    plt.title(f"PCA of ALT {args.encoding} ({X.shape[0]} sites)")
    plt.legend()
    plt.tight_layout()
    plt.savefig(f"{args.outprefix}.png", dpi=300)
    plt.savefig(f"{args.outprefix}.pdf")
    parser.report(sys.stderr)
    print(f"Sites used: {X.shape[0]}")
    print(f"Wrote: {args.outprefix}.scores.tsv, .png, .pdf")


if __name__ == "__main__":
    main()
