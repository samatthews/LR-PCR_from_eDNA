#!/bin/bash
# ---------------------------------------------------------------------------
# SLURM job array: one sample per job, each running collab_blast.sh.
#
# Edit the four lines marked EDIT, then:
#   1. pilot (200 reads, to measure speed):   see PILOT below
#   2. full run:   sbatch --array=1-$(ls fasta/*.fasta | wc -l) slurm_blast_array.sh
# Results: blast_out/<sample>.blast.tsv.gz, blast_out/<sample>.blast.time (time, RAM),
#          logs in slurm_logs/. Samples already finished are skipped if you resubmit.
#
# PILOT (run once before the full array):
#   mkdir -p pilot && awk '/^>/{n++} n<=200' fasta/barcode20.fasta > pilot/barcode20_pilot.fasta
#   sbatch --array=1 --time=04:00:00 --export=ALL,QUERY_DIR=pilot,OUT_DIR=pilot_out slurm_blast_array.sh
#   then:  grep -E "Elapsed|Maximum resident" pilot_out/barcode20_pilot.blast.time
# ---------------------------------------------------------------------------
#SBATCH --job-name=lrpcr_blast
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --output=slurm_logs/%x_%A_%a.out
##SBATCH --partition=EDIT_partition        # EDIT (remove one leading # to use)
##SBATCH --account=EDIT_account            # EDIT (remove one leading # to use)
set -euo pipefail

module load blast-plus 2>/dev/null || module load blast 2>/dev/null || true   # EDIT if your module is named differently
DB="${DB:-/path/to/core_nt/core_nt}"      # EDIT: database prefix (no extension)
QUERY_DIR="${QUERY_DIR:-fasta}"           # folder with <sample>.fasta from step 2
OUT_DIR="${OUT_DIR:-blast_out}"
SCRIPT_DIR="${SCRIPT_DIR:-$SLURM_SUBMIT_DIR}"   # folder holding collab_blast.sh

mkdir -p slurm_logs "$OUT_DIR"
FASTA=$(ls "$QUERY_DIR"/*.fasta | sort | sed -n "${SLURM_ARRAY_TASK_ID}p")
[ -n "$FASTA" ] || { echo "No sample number $SLURM_ARRAY_TASK_ID in $QUERY_DIR"; exit 1; }
S=$(basename "$FASTA" .fasta)
echo "Sample $S on $(hostname), $SLURM_CPUS_PER_TASK CPUs"

# collab_blast.sh works on a folder, so give it a folder holding just this sample.
ONE="$OUT_DIR/.query_$S"; mkdir -p "$ONE"
ln -sf "$(cd "$(dirname "$FASTA")" && pwd)/$(basename "$FASTA")" "$ONE/$S.fasta"
DB="$DB" THREADS="$SLURM_CPUS_PER_TASK" WORD_SIZE="${WORD_SIZE:-64}" \
  bash "$SCRIPT_DIR/collab_blast.sh" "$ONE" "$OUT_DIR"
rm -rf "$ONE"
