#!/bin/bash
# ---------------------------------------------------------------------------
# collab_blast.sh - BLAST step of the LR-PCR pipeline (step 03), standalone.
#
# Runs the same blastn search as 03_blast_reads.sh, on a different computer,
# and writes output that drops straight back into the pipeline (step 04/05).
#
# INPUT   One FASTA per sample, named <sample>.fasta (e.g. barcode20.fasta).
#         These are the step-02 outputs (reads/ folder): already Q- and
#         length-filtered. Do not rename the reads; the read IDs must match
#         the <sample>.readinfo.tsv files kept on the pipeline side.
#
# DATABASE  Any BLAST v5 nucleotide database with taxids (vert_mito from step 1,
#         or e.g. core_nt on a server). Use ONE database for every sample and
#         both assays, and record its name and date for the Methods. With a
#         database other than vert_mito, run 04a_lineages_from_blast.sh on the
#         returned files before steps 4/5.
#
# OUTPUT  One file per sample: <sample>.blast.tsv.gz
#         gzipped, tab-separated, no header, BLAST -outfmt 6 with these
#         16 columns in this order:
#           qseqid sseqid staxids pident length mismatch gapopen qstart qend
#           sstart send evalue bitscore qlen slen qcovs
#         Also writes blast_version.txt and compute_summary.tsv.
#         Send back the whole output folder; it goes into the pipeline's
#         blast/ folder (or subsample_N/blast/ for a subsample run).
#
# USAGE   DB=/path/to/refdb/vert_mito bash collab_blast.sh reads/ blast/
#         Optional: THREADS=32  WORD_SIZE=28
#           WORD_SIZE 64 = long-range PCR reads (8-10 kb), as in the paper
#           WORD_SIZE 28 = MiFish reads (~170 bp)
# ---------------------------------------------------------------------------
set -euo pipefail

QUERY_DIR="${1:?usage: DB=/path/to/vert_mito bash collab_blast.sh <query_dir> [out_dir]}"
OUT_DIR="${2:-blast}"
DB="${DB:?set DB to the vert_mito database prefix, e.g. DB=/data/refdb/vert_mito}"
THREADS="${THREADS:-8}"

# --- Search settings (identical to the pipeline's config.sh) ---------------
TASK=megablast
PERC_IDENTITY=85        # loose retrieval; species thresholds are applied later
WORD_SIZE="${WORD_SIZE:-64}"
EVALUE=1e-20
MAX_TARGET_SEQS=100     # enough subjects for the LCA to see congeners
FIELDS="qseqid sseqid staxids pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen qcovs"
# ---------------------------------------------------------------------------

[ -f "$DB.nal" ] || [ -f "$DB.nsq" ] || [ -f "$DB.00.nsq" ] || [ -f "$DB.njs" ] || {
  echo "No BLAST database found at $DB" >&2; exit 1; }
export BLASTDB="$(cd "$(dirname "$DB")" && pwd)"   # where taxdb.btd/.bti must sit
# staxids come from the database itself (BLAST db version 5, e.g. core_nt or vert_mito);
# the sanity check below warns if they come back empty.

mkdir -p "$OUT_DIR"
blastn -version > "$OUT_DIR/blast_version.txt"

# -mt_mode 1 (split by query) needs BLAST+ 2.12 or later; use it if available.
MT=(); if blastn -help 2>/dev/null | grep -q -- "-mt_mode"; then MT=(-mt_mode 1); fi

# Peak memory / wall time per sample: GNU time (Linux) -v, BSD time (macOS) -l.
if [ ! -x /usr/bin/time ]; then TIMECMD=()
elif /usr/bin/time -v true >/dev/null 2>&1; then TIMECMD=(/usr/bin/time -v -o)
else TIMECMD=(/usr/bin/time -l -o); fi

shopt -s nullglob
QUERIES=("$QUERY_DIR"/*.fasta)
[ ${#QUERIES[@]} -gt 0 ] || { echo "No .fasta files in $QUERY_DIR" >&2; exit 1; }

for q in "${QUERIES[@]}"; do
  s="$(basename "$q" .fasta)"; out="$OUT_DIR/$s.blast.tsv"
  [ -s "$q" ] || { echo "skip $s (empty)"; continue; }
  [ -s "$out.gz" ] && { echo "skip $s (already done)"; continue; }
  echo "BLAST $s: $(grep -c '>' "$q") reads"
  if [ ${#TIMECMD[@]} -gt 0 ]; then T=("${TIMECMD[@]}" "$OUT_DIR/$s.blast.time"); else T=(); fi
  ${T[@]+"${T[@]}"} blastn -task "$TASK" -query "$q" -db "$DB" \
    -perc_identity "$PERC_IDENTITY" -word_size "$WORD_SIZE" -evalue "$EVALUE" \
    -max_target_seqs "$MAX_TARGET_SEQS" -num_threads "$THREADS" ${MT[@]+"${MT[@]}"} \
    -outfmt "6 $FIELDS" -out "$out.part"
  mv "$out.part" "$out"

  # Sanity check: 16 columns, and taxids filled in (the assignment step needs them).
  if [ -s "$out" ]; then
    ncol=$(head -1 "$out" | awk -F'\t' '{print NF}')
    [ "$ncol" -eq 16 ] || echo "WARNING: $s has $ncol columns, expected 16" >&2
    bad=$(head -1000 "$out" | awk -F'\t' '$3=="N/A" || $3=="0" || $3==""' | wc -l)
    [ "$bad" -eq 0 ] || echo "WARNING: $s has $bad of the first 1000 hits without a taxid (check taxdb)" >&2
  fi
  gzip -f "$out"
done

# Compute summary (wall time, peak RAM per sample)
python3 - "$OUT_DIR" <<'PY'
import glob, os, re, sys
rows = []
for f in sorted(glob.glob(os.path.join(sys.argv[1], "*.blast.time"))):
    t = open(f).read()
    wall = re.search(r"Elapsed \(wall clock\) time.*?: ([\d:.]+)", t) or re.search(r"([\d.]+) real", t)
    gnu = re.search(r"Maximum resident set size \(kbytes\): (\d+)", t)
    bsd = re.search(r"(\d+)\s+maximum resident set size", t)
    gb = int(gnu.group(1)) / 1e6 if gnu else (int(bsd.group(1)) / 1e9 if bsd else float("nan"))
    rows.append(f"{os.path.basename(f)[:-11]}\t{wall.group(1) if wall else 'NA'}\t{gb:.2f}")
if rows:
    open(os.path.join(sys.argv[1], "compute_summary.tsv"), "w").write(
        "sample\twall_time\tpeak_RAM_GB\n" + "\n".join(rows) + "\n")
PY
echo "Done. Send back the folder: $OUT_DIR"
