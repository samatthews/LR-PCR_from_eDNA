#!/bin/bash
# Step 2 - length/quality filter each sample and record per-read length and mean Q.
#
# Input: FASTQ_DIR = the MinKNOW pass folder with barcodeNN/ subfolders (e.g. fastq_pass_raw).
#        Each chunk may exist as .fastq.gz, .fastq and _sanitized.fastq; ONE copy per chunk
#        is used (lib_fastq.sh), so no read is counted twice. Non-barcode folders are skipped.
# Output per sample: reads/combined_raw/<sample>.fastq.gz (all raw reads, one copy each)
#                    reads/<sample>.fasta (reads passing BLAST_MIN_Q and length)
#                    reads/<sample>.readinfo.tsv (read id, length, mean Q) - lets step 4
#                    test stricter Q cut-offs without re-running BLAST
# Also writes results/read_counts.tsv: raw vs filtered reads per sample (R1 L238-240),
#             results/raw_file_choice.tsv: which file was used for every chunk, and why.
# Time: minutes.
set -euo pipefail
source "$(dirname "$0")/config.sh"
mkdir -p "$READ_DIR" "$RESULT_DIR"

source "$(dirname "$0")/lib_fastq.sh"
command -v seqkit >/dev/null || { echo "seqkit not found - run: conda activate lrpcr" >&2; exit 1; }

# Samples: barcodeNN subfolders (MinKNOW layout) or one FASTQ file per sample.
# Other subfolders (unclassified, clean_concat_*, earlier outputs) are skipped.
declare -a NAMES=() SOURCES=()
if ls -d "$FASTQ_DIR"/*/ >/dev/null 2>&1; then
  for d in "$FASTQ_DIR"/*/; do
    n=$(basename "$d")
    case "$n" in
      barcode[0-9]*) NAMES+=("$n"); SOURCES+=("${d%/}") ;;
      *) echo "  skipping folder $n (not a barcode folder)" >&2 ;;
    esac
  done
else
  for f in "$FASTQ_DIR"/*.fastq "$FASTQ_DIR"/*.fastq.gz "$FASTQ_DIR"/*.fq "$FASTQ_DIR"/*.fq.gz; do
    [ -e "$f" ] || continue
    n=$(basename "$f"); n=${n%.gz}; n=${n%.fastq}; n=${n%.fq}
    NAMES+=("$n"); SOURCES+=("$f")
  done
fi
[ ${#NAMES[@]} -gt 0 ] || { echo "No barcode folders or FASTQ files found in $FASTQ_DIR" >&2; exit 1; }

RAW_DIR="$READ_DIR/combined_raw"; mkdir -p "$RAW_DIR"
CHOICE="$RESULT_DIR/raw_file_choice.tsv"
printf "sample\tchunk\tfile_used\treads\tbases\tdamaged\tplain_text_named_gz\tcopies_in_folder\tnote\n" > "$CHOICE"
printf "sample\tchunks\tdamaged_chunks\traw_reads\traw_bases\tduplicate_ids_removed\tpass_reads\tpass_bases\n" \
  > "$RESULT_DIR/read_counts.tsv"

for i in "${!NAMES[@]}"; do
  n=${NAMES[$i]}; src=${SOURCES[$i]}
  raw="$RAW_DIR/$n.fastq.gz"
  if [ -d "$src" ]; then
    # One file per chunk (see lib_fastq.sh), complete records only, then drop any read ID
    # that still appears twice (a safety net; the count is reported).
    picks=$(fq_pick "$src")
    [ -n "$picks" ] || { echo "  $n: no FASTQ files" >&2; continue; }
    printf "%s\n" "$picks" | awk -v s="$n" '{print s "\t" $0}' >> "$CHOICE"
    printf "%s\n" "$picks" | while IFS="$(printf '\t')" read -r key f r b dmg pl nc note; do
      if [ "$dmg" = 0 ]; then gzip -dcf "$f"; else fq_complete "$f"; fi
    done | seqkit rmdup -n -o "$raw" 2> "$RAW_DIR/$n.rmdup.log"
    chunks=$(printf "%s\n" "$picks" | wc -l | tr -d ' ')
    damaged=$(printf "%s\n" "$picks" | awk -F'\t' '$5 == 1' | wc -l | tr -d ' ')
  else
    fq_complete "$src" | seqkit rmdup -n -o "$raw" 2> "$RAW_DIR/$n.rmdup.log"
    chunks=1; damaged=$(fq_check "$src" | cut -f3)
  fi
  dups=$(grep -Eo '[0-9]+ duplicated records removed' "$RAW_DIR/$n.rmdup.log" | grep -Eo '^[0-9]+' || echo 0)
  rawstats=$(seqkit stats -T "$raw" | awk -F'\t' 'NR == 2 {print $4 "\t" $5}')

  # First filter: length window and mean quality (BLAST_MIN_Q, Q20).
  seqkit seq -m "$MIN_LEN" -M "$MAX_LEN" -Q "$BLAST_MIN_Q" "$raw" -o "$READ_DIR/$n.fastq.gz"
  seqkit fx2tab -n -i -l -q "$READ_DIR/$n.fastq.gz" | cut -f1-3 > "$READ_DIR/$n.readinfo.tsv"
  seqkit fq2fa "$READ_DIR/$n.fastq.gz" -o "$READ_DIR/$n.fasta"
  pass=$(awk -F'\t' '{r++; b+=$2} END {print r+0"\t"b+0}' "$READ_DIR/$n.readinfo.tsv")
  printf "%s\t%s\t%s\t%s\t%s\t%s\n" "$n" "$chunks" "$damaged" "$rawstats" "$dups" "$pass" >> "$RESULT_DIR/read_counts.tsv"
  echo "$n: $(cut -f1 <<< "$pass") of $(cut -f1 <<< "$rawstats") reads passed ($chunks chunks, $damaged damaged, $dups duplicate IDs removed)"
done
column -t "$RESULT_DIR/read_counts.tsv" 2>/dev/null || cat "$RESULT_DIR/read_counts.tsv"
echo "Per-chunk file choices: $CHOICE"
echo "Combined raw reads per sample (one copy of each read): $RAW_DIR"
