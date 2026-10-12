#!/bin/bash
# Step 2 (MiFish) - length-filter the MiFish reads and write them in the layout the
# rest of the pipeline expects. Settings come from config_mifish.sh.
#
# Input:  MIFISH_DIR, one file per sample (.fasta/.fa/.fna, .fastq/.fq, optionally .gz)
# Output: mifish/reads/<sample>.fasta          reads inside MIN_LEN-MAX_LEN (and >= Q20 for FASTQ)
#         mifish/reads/<sample>.readinfo.tsv   read id, length, mean Q ("NA" for FASTA input)
#         mifish/results/read_counts.tsv       raw vs kept reads, and how many were too short/long
# Time: seconds to minutes.
set -euo pipefail
source "${LRPCR_CONFIG:-$(dirname "$0")/config_mifish.sh}"
[ -d "$MIFISH_DIR" ] || { echo "MIFISH_DIR not found: $MIFISH_DIR (edit config_mifish.sh)" >&2; exit 1; }
mkdir -p "$READ_DIR" "$RESULT_DIR"

COUNTS="$RESULT_DIR/read_counts.tsv"
printf "sample\tinput\traw_reads\traw_bases\tpass_reads\tpass_bases\ttoo_short\ttoo_long\tlow_q\n" > "$COUNTS"

# id <TAB> length, read straight from a FASTA (handles wrapped sequences).
fa_lengths () {
  { gzip -dcf "$1"; } | awk '/^>/ {if (id != "") print id "\t" l; id = substr($1, 2); l = 0; next}
                              {l += length($0)}
                              END {if (id != "") print id "\t" l}'
}

shopt -s nullglob
n_files=0
for f in "$MIFISH_DIR"/*; do
  base=$(basename "$f")
  case "$base" in
    *.fasta|*.fa|*.fna|*.fasta.gz|*.fa.gz|*.fna.gz) kind=fasta ;;
    *.fastq|*.fq|*.fastq.gz|*.fq.gz)                 kind=fastq ;;
    *) continue ;;
  esac
  s=${base%.gz}; s=${s%.*}
  n_files=$((n_files + 1))

  if [ "$kind" = fasta ]; then
    # Length of every input read, then keep those inside the window.
    fa_lengths "$f" > "$READ_DIR/$s.all_lengths.tsv"
    seqkit seq -m "$MIN_LEN" -M "$MAX_LEN" "$f" -o "$READ_DIR/$s.fasta"
    fa_lengths "$READ_DIR/$s.fasta" | awk -F'\t' '{print $1 "\t" $2 "\tNA"}' > "$READ_DIR/$s.readinfo.tsv"
    low_q=NA
  else
    seqkit fx2tab -n -i -l "$f" | cut -f1,2 > "$READ_DIR/$s.all_lengths.tsv"
    seqkit seq -m "$MIN_LEN" -M "$MAX_LEN" "$f" -o "$READ_DIR/$s.len.fastq.gz"
    seqkit seq -Q "$BLAST_MIN_Q" "$READ_DIR/$s.len.fastq.gz" -o "$READ_DIR/$s.fastq.gz"
    seqkit fx2tab -n -i -l -q "$READ_DIR/$s.fastq.gz" | cut -f1-3 > "$READ_DIR/$s.readinfo.tsv"
    seqkit fq2fa "$READ_DIR/$s.fastq.gz" -o "$READ_DIR/$s.fasta"
    in_len=$(seqkit stats -T "$READ_DIR/$s.len.fastq.gz" | awk -F'\t' 'NR == 2 {print $4}')
    low_q=$(( in_len - $(wc -l < "$READ_DIR/$s.readinfo.tsv") ))
    rm -f "$READ_DIR/$s.len.fastq.gz"
  fi

  # Read IDs must be unique within a sample (the assignment step joins on them).
  dup=$(cut -f1 "$READ_DIR/$s.readinfo.tsv" | sort | uniq -d | wc -l | tr -d ' ')
  [ "$dup" -eq 0 ] || echo "WARNING: $s has $dup duplicated read IDs" >&2

  stats=$(awk -F'\t' -v lo="$MIN_LEN" -v hi="$MAX_LEN" \
    '{r++; b += $2; if ($2 < lo) s++; else if ($2 > hi) l++} END {print r+0 "\t" b+0 "\t" s+0 "\t" l+0}' \
    "$READ_DIR/$s.all_lengths.tsv")
  read -r raw_r raw_b short long <<< "$stats"
  pass=$(awk -F'\t' '{r++; b += $2} END {print r+0 "\t" b+0}' "$READ_DIR/$s.readinfo.tsv")
  printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" "$s" "$kind" "$raw_r" "$raw_b" "$pass" "$short" "$long" "$low_q" >> "$COUNTS"
  rm -f "$READ_DIR/$s.all_lengths.tsv"
  echo "$s ($kind): kept $(cut -f1 <<< "$pass") of $raw_r reads ($short too short, $long too long)"
done
[ "$n_files" -gt 0 ] || { echo "No FASTA/FASTQ files in $MIFISH_DIR" >&2; exit 1; }
column -t "$COUNTS" 2>/dev/null || cat "$COUNTS"
