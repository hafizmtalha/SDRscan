#!/usr/bin/env python3
"""
Write a findZX units.tsv from the pipeline sample sheet.

findZX compares a "homogametic" and a "heterogametic" group. Which sex is
which depends on the system being tested:

  --system XY   males heterogametic, females homogametic
  --system ZW   females heterogametic, males homogametic

Using findZX's own group labels avoids editing workflow/rules/common.smk.

--shuffle produces a negative control: the same samples with sex labels
randomly reassigned (keeping the number of each sex), so any "sex-linked"
signal that survives is not driven by sex.
"""

import argparse
import random

import pandas as pd


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--samples", default="config/samples.tsv")
    ap.add_argument("--system", choices=["XY", "ZW"], required=True)
    ap.add_argument("--out", default="units.tsv")
    ap.add_argument("--shuffle", action="store_true", help="Randomly reassign sex labels (negative control)")
    ap.add_argument("--seed", type=int, default=1)
    args = ap.parse_args()

    df = pd.read_csv(args.samples, sep="\t")
    sex = df["sex"].tolist()
    if args.shuffle:
        random.Random(args.seed).shuffle(sex)
        df["sex"] = sex

    hetero = "M" if args.system == "XY" else "F"
    df["group"] = ["heterogametic" if s == hetero else "homogametic" for s in df["sex"]]
    df[["sample", "group", "fq1", "fq2"]].to_csv(args.out, sep="\t", index=False)

    print(df[["sample", "sex", "group"]].to_string(index=False))
    print(f"Wrote: {args.out}")


if __name__ == "__main__":
    main()
