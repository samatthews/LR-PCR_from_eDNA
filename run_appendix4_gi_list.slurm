#!/bin/bash
#SBATCH --job-name=lrpcr_gilist
#SBATCH --account=coenv
#SBATCH --partition=cpu-g2
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=01:00:00
#SBATCH --output=lrpcr_gilist_%j.out        # %j becomes the job ID
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=EDIT_YOUR_EMAIL

# Appendix 4 - every record in the vert_mtDNA BLAST database: GI, accession.version,
# taxid and sequence length, one row per record (about 2.2 million rows).
#   sbatch run_appendix4_gi_list.slurm
#
# Output ($OUT):
#   Appendix4_vert_mtDNA_records.tsv.gz   gi, accession, taxid, length (header row first)
#   Appendix4_vert_mtDNA_gi_list.txt.gz   GI numbers only, one per line
#   vert_mtDNA_alias_info.txt             the database's alias file and blastdbcmd -info,
#                                         to check against the record count

set -euo pipefail

SIF=EDIT_PATH_TO_SIF
DB=EDIT_PATH_TO_DB
OUT=EDIT_PATH_TO_OUT

mkdir -p "$OUT"

run () {
  apptainer exec --bind /gscratch/ --env BLASTDB="$(dirname "$DB")" "$SIF" "$@"
}

# Record of how the database was built: the alias file (.nal) names the GI list it was made
# from, and -info gives the number of sequences and the build date.
{
  echo "== alias file(s)"
  for f in "$DB".nal "$DB".pal; do [ -f "$f" ] && { echo "-- $f"; cat "$f"; }; done
  echo; echo "== blastdbcmd -info"
  run blastdbcmd -db "$DB" -info
} > "$OUT/vert_mtDNA_alias_info.txt" 2>&1 || true

echo "[$(date)] Listing records in $DB"
{
  printf "gi\taccession\ttaxid\tlength\n"
  run blastdbcmd -db "$DB" -entry all -outfmt $'%g\t%a\t%T\t%l'
} | gzip > "$OUT/Appendix4_vert_mtDNA_records.tsv.gz"

zcat "$OUT/Appendix4_vert_mtDNA_records.tsv.gz" | tail -n +2 | cut -f1 | gzip > "$OUT/Appendix4_vert_mtDNA_gi_list.txt.gz"

n=$(zcat "$OUT/Appendix4_vert_mtDNA_records.tsv.gz" | tail -n +2 | wc -l)
na=$(zcat "$OUT/Appendix4_vert_mtDNA_gi_list.txt.gz" | grep -c -v '^[0-9][0-9]*$' || true)
echo "[$(date)] $n records written; $na without a GI number (shown as N/A or 0)."
echo "Outputs in $OUT"
