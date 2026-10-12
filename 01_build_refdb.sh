#!/bin/bash
# Step 1 - build one local BLAST database from the 12S and COI GI lists.
#
# The original analysis aliased NCBI core_nt with these lists. Downloading core_nt
# is >150 GB, so instead we fetch just the listed records (~710k sequences, union of
# both lists) from NCBI and build a small database with taxonomy attached.
# One combined database is used for both fragments, so every read is compared with
# full mitogenomes AND partial 12S/COI records (R1 major comment 2).
#
# Time: ~1-3 h, mostly waiting on NCBI (faster with NCBI_API_KEY). Safe to re-run:
# finished chunks are skipped, so an interrupted download resumes where it stopped.
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"   # resolve before any cd
source "$SCRIPT_DIR/config.sh"

# Check tools first (most errors here mean the conda env is not active).
for t in seqkit makeblastdb taxonkit python3 curl; do
  command -v "$t" >/dev/null || { echo "'$t' not found. Run: conda activate lrpcr   (or bash 00_setup_env.sh)" >&2; exit 1; }
done

mkdir -p "$REF_DIR/chunks" "$LOG_DIR"
cd "$REF_DIR"
rm -f chunks/gi_* chunks/.split_done   # leftovers from the earlier EDirect version

# 1. Union of the two lists --------------------------------------------------
cat "$GILIST_12S" "$GILIST_COI" | tr -d '\r' | grep -E '^[0-9]+$' | sort -u > all_gis.txt
echo "$(wc -l < all_gis.txt) unique GIs in the two lists"

# 2. Fetch sequences + metadata from NCBI (resumable; no EDirect needed) -----
fetch_status=0
python3 "$SCRIPT_DIR/fetch_ncbi.py" all_gis.txt chunks || fetch_status=$?
if [ "$fetch_status" -ne 0 ]; then
  if [ "${ALLOW_PARTIAL:-0}" = 1 ] && ls chunks/c_*.tsv >/dev/null 2>&1; then
    echo "ALLOW_PARTIAL=1: building the database from the batches downloaded so far." >&2
    echo "  Missing GIs are in refdb/chunks/missing_gis.txt. Re-run without ALLOW_PARTIAL later to" >&2
    echo "  fill the gaps; the database is rebuilt from all batches each time." >&2
  else
    echo >&2
    echo "Some batches did not download. Either re-run this script (only missing batches are fetched)," >&2
    echo "or build now from what you have:   ALLOW_PARTIAL=1 bash 01_build_refdb.sh" >&2
    exit 1
  fi
fi
if [ -s chunks/unavailable_records.tsv ]; then
  echo "$(cut -f2 chunks/unavailable_records.tsv | sort -u | wc -l) records listed by NCBI but not served (removed/replaced);"
  echo "  skipped and listed in refdb/chunks/unavailable_records.tsv"
fi

# 3. Combine, drop duplicate accessions -------------------------------------
cat chunks/c_*.fa | seqkit rmdup -n -o refs.fasta 2> "$LOG_DIR/rmdup.log"
cat chunks/c_*.tsv | awk -F'\t' 'NF>=2 && !seen[$1]++' > ref_metadata.tsv
cut -f1,2 ref_metadata.tsv > acc2taxid.txt
seqkit stats -T refs.fasta | tee "$LOG_DIR/refdb_stats.tsv"

# 4. BLAST database with taxids ---------------------------------------------
rm -f vert_mito.*
makeblastdb -in refs.fasta -dbtype nucl -parse_seqids -blastdb_version 5 \
  -taxid_map acc2taxid.txt -out vert_mito \
  -title "Vertebrate mito 12S+COI+mitogenomes (from core_nt GI lists)" \
  | tee "$LOG_DIR/makeblastdb.log"

# 5. Taxonomy names for BLAST (taxdb) and lineages for the LCA step ------------
if [ ! -f taxdb.btd ]; then
  curl -sSLO https://ftp.ncbi.nlm.nih.gov/blast/db/taxdb.tar.gz && tar xzf taxdb.tar.gz && rm taxdb.tar.gz
fi
TAXDIR="$REF_DIR/taxdump"
if [ ! -f "$TAXDIR/nodes.dmp" ]; then
  mkdir -p "$TAXDIR"
  curl -sSL https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz | tar xz -C "$TAXDIR"
fi
cut -f2 acc2taxid.txt | sort -u |
  taxonkit reformat --data-dir "$TAXDIR" -I 1 -f "{k};{p};{c};{o};{f};{g};{s}" > taxid_lineage.tsv

echo
echo "Database ready: $REF_DIR/vert_mito  ($(wc -l < ref_metadata.tsv) records, $(wc -l < taxid_lineage.tsv) taxa)"
