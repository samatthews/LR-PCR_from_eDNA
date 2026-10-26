#!/bin/bash
# Step 4a - only needed when BLAST was run against a database OTHER than refdb/vert_mito
# (e.g. core_nt on the lab server, via collab_blast.sh).
#
# The assignment step (lrpcr_assign.py) looks up each hit's taxid in
# refdb/taxid_lineage.tsv. That file was built for the vert_mito references only, so
# taxids that exist only in the other database would have no lineage, and their reads
# would end up unassigned. This script adds a lineage for every taxid found in the BLAST
# results, using the NCBI taxonomy in refdb/taxdump (downloaded by step 1).
#
# Usage:  bash 04a_lineages_from_blast.sh            (reads every file in BLAST_DIR)
#         bash 04a_lineages_from_blast.sh other_dir  (a different folder of .blast.tsv.gz)
# Run it again whenever new BLAST results are added. The original file is kept as
# taxid_lineage.vert_mito.tsv the first time.
set -euo pipefail
source "${LRPCR_CONFIG:-$(dirname "$0")/config.sh}"
SRC="${1:-$BLAST_DIR}"
TAXDIR="$REF_DIR/taxdump"
LIN="$REF_DIR/taxid_lineage.tsv"

command -v taxonkit >/dev/null || { echo "taxonkit not found - run: conda activate lrpcr" >&2; exit 1; }
[ -f "$TAXDIR/nodes.dmp" ] || { echo "No NCBI taxonomy in $TAXDIR - run 01_build_refdb.sh, or download" \
  "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz and unpack it there" >&2; exit 1; }
ls "$SRC"/*.blast.tsv* >/dev/null 2>&1 || { echo "No .blast.tsv(.gz) files in $SRC" >&2; exit 1; }

[ -f "$LIN" ] && [ ! -f "$REF_DIR/taxid_lineage.vert_mito.tsv" ] && cp "$LIN" "$REF_DIR/taxid_lineage.vert_mito.tsv"
touch "$LIN"

# Every taxid in column 3 (staxids; several taxids are separated by ';' - all are kept).
TMP="$REF_DIR/.taxids_from_blast.txt"
for f in "$SRC"/*.blast.tsv*; do gzip -dcf "$f"; done | cut -f3 | tr ';' '\n' |
  grep -E '^[0-9]+$' | grep -vx 0 | sort -u > "$TMP"
NEW="$REF_DIR/.taxids_new.txt"
cut -f1 "$LIN" | sort -u | comm -23 "$TMP" - > "$NEW"
echo "$(wc -l < "$TMP" | tr -d ' ') taxids in the BLAST results; $(wc -l < "$NEW" | tr -d ' ') not yet in taxid_lineage.tsv"

if [ -s "$NEW" ]; then
  taxonkit reformat --data-dir "$TAXDIR" -I 1 -f "{k};{p};{c};{o};{f};{g};{s}" < "$NEW" > "$NEW.lin"
  # taxids missing from this taxdump (newer than the download) come back with an empty lineage
  miss=$(awk -F'\t' '$2 == "" || $2 == ";;;;;;"' "$NEW.lin" | wc -l | tr -d ' ')
  cat "$NEW.lin" >> "$LIN"
  echo "Added $(wc -l < "$NEW.lin" | tr -d ' ') lineages to $LIN"
  if [ "$miss" -gt 0 ]; then
    awk -F'\t' '$2 == "" || $2 == ";;;;;;" {print $1}' "$NEW.lin" > "$REF_DIR/taxids_without_lineage.txt"
    echo "WARNING: $miss taxids are not in $TAXDIR (taxonomy older than the database?)." >&2
    echo "  Listed in $REF_DIR/taxids_without_lineage.txt. To fix: delete $TAXDIR, download a" >&2
    echo "  fresh taxdump.tar.gz into it, delete $LIN, restore taxid_lineage.vert_mito.tsv as" >&2
    echo "  taxid_lineage.tsv, and run this script again." >&2
  fi
  rm -f "$NEW.lin"
fi
rm -f "$TMP" "$NEW"
