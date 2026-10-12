#!/bin/bash
# Count reads and bases per barcode in the ORIGINAL MinKNOW output (pass and fail),
# for R1 L238-240 (total reads generated / discarded) and 08_read_accounting.py.
#
# Each chunk is counted ONCE even when the folder holds several copies of it
# (.fastq.gz, .fastq, _sanitized.fastq); see lib_fastq.sh for which copy is used.
# A file cut off at the end of the run is counted up to the damage.
# Only barcodeNN/ and unclassified/ subfolders are counted (clean_concat_* etc. are skipped).
#
# Usage:  bash count_raw_reads.sh <run folder> > raw_read_totals.tsv
#   e.g.  bash count_raw_reads.sh ~/.../20251120_2032_MD-100902_FBE54526_4090b961 > raw_read_totals.tsv
# The run folder is the one containing fastq_pass*/ (e.g. fastq_pass_raw/) and fastq_fail*/.
set -uo pipefail
RUN="${1:?give the run folder that contains fastq_pass*/ and fastq_fail*/}"
source "$(dirname "$0")/lib_fastq.sh"

# Two pass folders (e.g. fastq_pass and fastq_pass_raw) would usually be copies of each other.
n_pass=$(ls -d "$RUN"/fastq_pass*/ 2>/dev/null | wc -l | tr -d ' ')
[ "$n_pass" -le 1 ] || echo "WARNING: $n_pass fastq_pass* folders in $RUN - all are counted; if they are copies, move the extra one out first" >&2

printf "folder\tbarcode\treads\tbases\tfiles\tplain_text_named_gz\ttruncated_files\tduplicate_copies_ignored\n"
for d in "$RUN"/fastq_pass*/*/ "$RUN"/fastq_fail*/*/; do
  [ -d "$d" ] || continue
  bc=$(basename "$d")
  case "$bc" in barcode[0-9]*|unclassified) ;; *) echo "  skipping $(basename "$(dirname "$d")")/$bc" >&2; continue ;; esac
  fq_pick "${d%/}" | awk -F'\t' -v folder="$(basename "$(dirname "$d")")" -v bc="$bc" '
    { r += $3; b += $4; n++; dmg += $5; pl += $6; extra += $7 - 1 }
    END { printf "%s\t%s\t%d\t%d\t%d\t%d\t%d\t%d\n", folder, bc, r, b, n, pl, dmg, extra }'
done
