#!/usr/bin/env python3
"""
PCA on ALT-presence genotypes in a region, coloured by recorded sex.

Each genotype becomes 1 (ALT present), 0 (REF only) or missing. Sites with
fewer than --min-called called samples are dropped; remaining missing values
are imputed with the site mean. Samples are labelled but sex labels are not
used to fit the PCA, so a sample clustering with the other sex is a hint of
mis-sexing or of no sex-linked structure in the region.
"""

import argparse
import re

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA

MISSING = {".", "./.", ".|.", "./././.", ".|.|.|.", "././././.", ".|.|.|.|."}


def read_list(path):
    with open(path) as f:
        return [x.strip() for x in f if x.strip()]


def alt_presence(gt):
    if gt is None or gt.strip() in MISSING:
        return -1
    return int(any(a not in {"", ".", "0"} for a in re.split(r"[\/|]", gt.strip())))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--genotypes", required=True)
    ap.add_argument("--vcf-samples", default="vcf_samples.txt")
    ap.add_argument("--males", default="males.txt")
    ap.add_argument("--females", default="females.txt")
    ap.add_argument("--outprefix", default="region_PCA")
    ap.add_argument("--min-called", type=int, default=None,
                    help="Min called samples per site (default: 80%% of samples)")
    ap.add_argument("--min-sites", type=int, default=100)
    ap.add_argument("--max-sites", type=int, default=None,
                    help="Randomly subsample to at most this many sites")
    ap.add_argument("--seed", type=int, default=1)
    args = ap.parse_args()

    samples = read_list(args.vcf_samples)
    males = set(read_list(args.males))
    females = set(read_list(args.females))
    min_called = args.min_called if args.min_called is not None else int(0.8 * len(samples))

    rows = []
    with open(args.genotypes) as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            vals = [alt_presence(gt) for gt in parts[2:]]
            called = [v for v in vals if v >= 0]
            if len(called) < min_called:
                continue
            mean = float(np.mean(called))
            rows.append([mean if v < 0 else v for v in vals])

    X = np.array(rows, dtype=float)
    if args.max_sites and X.shape[0] > args.max_sites:
        idx = np.random.default_rng(args.seed).choice(X.shape[0], args.max_sites, replace=False)
        X = X[np.sort(idx)]
    if X.shape[0] < args.min_sites:
        raise SystemExit(f"Too few sites after filtering: {X.shape[0]}")

    X = X - X.mean(axis=1, keepdims=True)
    pca = PCA(n_components=2)
    pcs = pca.fit_transform(X.T)

    out = pd.DataFrame({
        "sample": samples,
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
    plt.title(f"PCA of ALT presence ({X.shape[0]} sites)")
    plt.legend()
    plt.tight_layout()
    plt.savefig(f"{args.outprefix}.png", dpi=300)
    plt.savefig(f"{args.outprefix}.pdf")
    print(f"Sites used: {X.shape[0]}")
    print(f"Wrote: {args.outprefix}.scores.tsv, .png, .pdf")


if __name__ == "__main__":
    main()
