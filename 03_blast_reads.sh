#!/bin/bash
# Step 3 - BLAST every filtered read against the combined reference database.
# Adapted from longread_blast20251211_12s_exp.sh, with settings for raw reads:
#   - loose retrieval (PERC_IDENTITY 85, WORD_SIZE 64): the species threshold is set
#     later from the data, not by BLAST
#   - 100 target sequences instead of 5, so the LCA sees congeners
#   - taxids and subject length reported for the assignment step
# Runtime and peak memory are logged per sample to logs/ (answers R1 L81-84).
# Usage: bash 03_blast_reads.sh              (all samples)
#        bash 03_blast_reads.sh --aquarium   (the AQ_SAMPLES in config.sh)
#        bash 03_blast_reads.sh --aquarium --subsample 5000   (quick calibration pilot)
#        bash 03_blast_reads.sh barcode19     (one or more named samples)
#        LRPCR_CONFIG=config_mifish.sh bash 03_blast_reads.sh   (MiFish reads)
set -euo pipefail
# LRPCR_CONFIG selects another config, e.g. LRPCR_CONFIG=config_mifish.sh for the MiFish reads.
source "${LRPCR_CONFIG:-$(dirname "$0")/config.sh}"
mkdir -p "$BLAST_DIR" "$LOG_DIR"

DB="$REF_DIR/vert_mito"
[ -f "$DB.njs" ] || [ -f "$DB.nal" ] || { echo "No database at $DB - run 01_build_refdb.sh" >&2; exit 1; }

# Options:
#   --aquarium        the AQ_SAMPLES list from config.sh
#   --subsample N     BLAST a random N reads per sample (seed 11) into subsample_N/;
#                     used for the quick calibration (05_calibrate.py --subsample N)
# Named samples are split on spaces, so a list passed as one quoted word (zsh does this
# with $AQ_SAMPLES) still works.
SUB=""; ARGS=()
while [ $# -gt 0 ]; do
  case "$1" in
    --aquarium) [ -n "$AQ_SAMPLES" ] || { echo "AQ_SAMPLES is empty in config.sh" >&2; exit 1; }
                for s in $AQ_SAMPLES; do ARGS+=("$s"); done ;;
    --subsample) SUB="${2:?--subsample needs a number}"; shift ;;
    *) for s in $1; do ARGS+=("$s"); done ;;
  esac
  shift
done
SAMPLES=()
if [ ${#ARGS[@]} -gt 0 ]; then SAMPLES=("${ARGS[@]}"); else
  for f in "$READ_DIR"/*.fasta; do SAMPLES+=("$(basename "$f" .fasta)"); done
fi

# Subsample mode: make a small reads folder (FASTA + matching readinfo) per sample.
if [ -n "$SUB" ]; then
  SUBDIR="$(dirname "$READ_DIR")/subsample_$SUB"   # beside reads/: PROJECT_DIR for LR-PCR, mifish/ for MiFish
  mkdir -p "$SUBDIR/reads" "$SUBDIR/blast"
  for s in "${SAMPLES[@]}"; do
    src="$READ_DIR/$s.fasta"; dst="$SUBDIR/reads/$s.fasta"
    [ -s "$src" ] || continue
    if [ ! -s "$SUBDIR/reads/$s.readinfo.tsv" ]; then
      # Shuffle to a file first: piping into head cuts the shuffle off, which the strict
      # error mode (pipefail) treats as a failure and silently stops the script.
      seqkit shuffle -s 11 "$src" -o "$dst.shuf" 2>/dev/null
      seqkit head -n "$SUB" "$dst.shuf" -o "$dst"
      rm -f "$dst.shuf"
      grep '>' "$dst" | sed 's/^>//; s/[[:space:]].*//' > "$dst.ids"
      awk -F'\t' 'NR==FNR {keep[$1]=1; next} ($1 in keep)' "$dst.ids" "$READ_DIR/$s.readinfo.tsv" \
        > "$SUBDIR/reads/$s.readinfo.tsv"
      rm -f "$dst.ids"
      echo "  $s: $(grep -c '>' "$dst") reads sampled"
    fi
  done
  READ_DIR="$SUBDIR/reads"; BLAST_DIR="$SUBDIR/blast"
  echo "Subsample mode: $SUB random reads per sample -> $SUBDIR"
fi

# GNU time (Linux) uses -v; BSD time (macOS) uses -l. Both report peak memory.
if [ ! -x /usr/bin/time ]; then TIMECMD=()
elif /usr/bin/time -v true >/dev/null 2>&1; then TIMECMD=(/usr/bin/time -v -o)
else TIMECMD=(/usr/bin/time -l -o); fi

FIELDS="qseqid sseqid staxids pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen qcovs"
for s in "${SAMPLES[@]}"; do
  q="$READ_DIR/$s.fasta"; out="$BLAST_DIR/$s.blast.tsv"
  [ -s "$q" ] || { echo "skip $s (no reads file at $q)"; continue; }
  [ -s "$out.gz" ] && { echo "skip $s (already done)"; continue; }
  echo "BLAST $s: $(grep -c '>' "$q") reads"
  if [ ${#TIMECMD[@]} -gt 0 ]; then T=("${TIMECMD[@]}" "$LOG_DIR/$s.blast.time"); else T=(); fi
  ${T[@]+"${T[@]}"} blastn -task megablast -query "$q" -db "$DB" \
      -perc_identity "$PERC_IDENTITY" -word_size "$WORD_SIZE" -evalue "$EVALUE" \
      -max_target_seqs "$MAX_TARGET_SEQS" -num_threads "$THREADS" -mt_mode 1 \
      -outfmt "6 $FIELDS" -out "$out.part"
  mv "$out.part" "$out"
  gzip -f "$out"
done

# One-line compute summary per sample (wall time, peak RAM) for the Methods.
python3 - "$LOG_DIR" <<'PY'
import glob, os, re, sys
rows = []
for f in sorted(glob.glob(os.path.join(sys.argv[1], "*.blast.time"))):
    t = open(f).read()
    wall = re.search(r"Elapsed \(wall clock\) time.*?: ([\d:.]+)", t) or re.search(r"([\d.]+) real", t)
    mem_gnu = re.search(r"Maximum resident set size \(kbytes\): (\d+)", t)
    mem_bsd = re.search(r"(\d+)\s+maximum resident set size", t)
    gb = int(mem_gnu.group(1)) / 1e6 if mem_gnu else (int(mem_bsd.group(1)) / 1e9 if mem_bsd else float("nan"))
    rows.append(f"{os.path.basename(f).replace('.blast.time','')}\t{wall.group(1) if wall else 'NA'}\t{gb:.2f}")
if rows:
    out = os.path.join(sys.argv[1], "compute_summary.tsv")
    open(out, "w").write("sample\twall_time\tpeak_RAM_GB\n" + "\n".join(rows) + "\n")
    print("\n".join(["sample\twall_time\tpeak_RAM_GB"] + rows))
PY
