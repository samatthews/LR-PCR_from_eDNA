#!/bin/bash
#SBATCH --job-name=lrpcr_refcov
#SBATCH --account=coenv
#SBATCH --partition=cpu-g2
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=lrpcr_refcov_%j.out        # %j becomes the job ID
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=EDIT_YOUR_EMAIL

# Step 6 - reference coverage of the aquarium species in vert_mtDNA.
# Submit from the lrpcr_blast folder (it needs 06_reference_coverage.py and
# lists/aquarium_species.tsv):   sbatch run_06_reference_coverage.slurm
#
# 1. blastdbcmd dumps one line per database record: accession, taxid, length,
#    scientific name, title  ->  $OUT/ref_metadata.tsv
# 2. 06_reference_coverage.py counts complete mitogenomes, other records >= 1 kb,
#    12S and COI records for each tank species  ->  $OUT/reference_coverage.tsv

set -euo pipefail

SIF=EDIT_PATH_TO_SIF
DB=EDIT_PATH_TO_DB
OUT=EDIT_PATH_TO_OUT

cd "${SLURM_SUBMIT_DIR:-$PWD}"
mkdir -p "$OUT"

# BLASTDB points blastdbcmd at the database folder, where the BLAST taxonomy files
# (taxdb.btd/.bti) live if they were installed; without them the scientific name is
# "N/A" and the script takes the species from the record title instead.
run () {
  apptainer exec --bind /gscratch/ --env BLASTDB="$(dirname "$DB")" "$SIF" "$@"
}

echo "[$(date)] Dumping record metadata from $DB"
run blastdbcmd -db "$DB" -entry all -outfmt $'%a\t%T\t%l\t%S\t%t' > "$OUT/ref_metadata.tsv"
echo "[$(date)] $(wc -l < "$OUT/ref_metadata.tsv") records written to $OUT/ref_metadata.tsv"

if run python3 -c "import pandas" 2>/dev/null; then
  run python3 06_reference_coverage.py \
    --metadata "$OUT/ref_metadata.tsv" \
    --species lists/aquarium_species.tsv \
    --out "$OUT/reference_coverage.tsv"
  echo "[$(date)] Done: $OUT/reference_coverage.tsv"
else
  echo "The container has no Python with pandas. Copy $OUT/ref_metadata.tsv to your computer and run:"
  echo "  python3 06_reference_coverage.py --metadata ref_metadata.tsv --out reference_coverage.tsv"
fi
