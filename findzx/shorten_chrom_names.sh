#!/usr/bin/env bash
# Shorten long chromosome names in findZX output tables before replotting,
# e.g. "Species_name_female_chr_07" -> "chr_07", and sort rows by chromosome.
#
# Usage (inside results/<run>/output/no_synteny/tables or a copy of it):
#   findzx/shorten_chrom_names.sh 'Species_name_female_chr_' 'chr_' file1 file2 ...
#
# Then regenerate only the plots:
#   snakemake -s workflow/findZX --configfile <cfg> --cores 4 --rerun-triggers mtime \
#     -- results/<run>/output/no_synteny/plots/.misc/plotting_chr.done
set -euo pipefail

FROM="${1:?usage: $0 FROM_PREFIX TO_PREFIX files...}"
TO="$2"
shift 2

for f in "$@"; do
  sed -i -E "s/${FROM}/${TO}/g" "$f"
  # keep a header line (if the first line has no numeric 2nd column) on top
  if head -n 1 "$f" | awk -F'\t' '{exit ($2 ~ /^[0-9]+$/)}'; then
    { head -n 1 "$f"; tail -n +2 "$f" | sort -k1,1V -k2,2n; } > "$f.tmp"
  else
    sort -k1,1V -k2,2n "$f" > "$f.tmp"
  fi
  mv "$f.tmp" "$f"
done
