"""
Shared genotype helpers for any ploidy.

Genotypes come from `bcftools query -f '%CHROM\\t%POS[\\t%GT]\\n'` and look like
0/1, 0/0/0/1, 0|1|1|1, 1/1/1/1/1/1, ./. etc. Ploidy is declared by the user
(PLOIDY in the pipeline config, optionally per sample) rather than guessed, and
genotypes whose allele count does not match the declared ploidy are treated as
missing and counted, so a mismatch between calling and analysis is visible.
"""

import math
import re

_SPLIT = re.compile(r"[\/|]")


def read_list(path):
    with open(path) as f:
        return [x.strip() for x in f if x.strip()]


def load_ploidy(samples, ploidy=None, sample_ploidy_file=None):
    """Return {sample: ploidy}. A per-sample file (sample<TAB>ploidy) overrides --ploidy."""
    table = {}
    if sample_ploidy_file:
        with open(sample_ploidy_file) as f:
            for line in f:
                parts = line.split()
                if len(parts) >= 2 and parts[1].isdigit():
                    table[parts[0]] = int(parts[1])
    out = {}
    for s in samples:
        p = table.get(s, ploidy)
        if p is None:
            raise SystemExit(f"ERROR: no ploidy for sample {s}; give --ploidy or --sample-ploidy")
        if p < 1:
            raise SystemExit(f"ERROR: invalid ploidy {p} for sample {s}")
        out[s] = int(p)
    return out


class GenotypeParser:
    """Turns GT strings into ALT dosage (number of non-reference alleles) or None if missing."""

    def __init__(self, ploidy_by_sample=None):
        self.ploidy = ploidy_by_sample or {}
        self.mismatched = 0
        self.mismatch_examples = {}

    def dosage(self, gt, sample=None):
        if gt is None:
            return None
        alleles = _SPLIT.split(gt.strip())
        if all(a in {"", "."} for a in alleles):
            return None
        expected = self.ploidy.get(sample)
        if expected is not None and len(alleles) != expected:
            self.mismatched += 1
            self.mismatch_examples.setdefault(sample, gt)
            return None
        if any(a in {"", "."} for a in alleles):
            return None  # partially missing call
        return sum(1 for a in alleles if a != "0")

    def report(self, stream):
        if self.mismatched:
            ex = ", ".join(f"{s}: {g}" for s, g in list(self.mismatch_examples.items())[:5])
            print(f"WARNING: {self.mismatched} genotypes did not match the declared ploidy and were "
                  f"treated as missing (e.g. {ex}). Check PLOIDY / the ploidy column.", file=stream)


def threshold(value, n):
    """A threshold given as a count (7) or a fraction of n samples (0.8 -> ceil(0.8*n))."""
    v = float(value)
    if v < 0:
        raise SystemExit(f"ERROR: negative threshold {value}")
    if "." in str(value) and v <= 1:
        return int(math.ceil(v * n - 1e-9))
    return int(v)


def check_dosages(dosages, ploidy_by_sample):
    if not dosages:
        return
    pmax = max(ploidy_by_sample.values())
    bad = [d for d in dosages if d < 1 or d > pmax]
    if bad:
        raise SystemExit(f"ERROR: --alt-dosage {bad} outside 1..{pmax} (max ploidy)")
