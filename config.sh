#!/bin/bash
# ---------------------------------------------------------------------------
# LR-PCR BLAST-only pipeline: settings shared by every step.
# Edit the paths in section 1, then run the numbered scripts in order.
# ---------------------------------------------------------------------------

# 1. PATHS -------------------------------------------------------------------
PROJECT_DIR="/Users/stephaniematthews/Documents/Rockfish/MinION_steph/MolEcolRevisions/lrpcr_blast"                 # this folder
FASTQ_DIR="/Users/stephaniematthews/Documents/Rockfish/MinION_steph/LRPCR_Manuscript_Analysis/LRPCR_20251120/no_sample_id/20251120_2032_MD-100902_FBE54526_4090b961/fastq_pass_raw"             # MinKNOW pass folder with barcodeNN/ subfolders
                                                # (duplicate copies of a chunk are used once;
                                                # non-barcode folders are skipped), or one
                                                # .fastq/.fastq.gz file per sample
GILIST_12S="$PROJECT_DIR/lists/vert_12s_target_gilist.txt"
GILIST_COI="$PROJECT_DIR/lists/vert_COI_target_gilist.txt"

REF_DIR="$PROJECT_DIR/refdb"                    # BLAST database is built here (~1-2 GB)
READ_DIR="$PROJECT_DIR/reads"                   # filtered reads (FASTA) + per-read stats
BLAST_DIR="$PROJECT_DIR/blast"                  # raw BLAST tables (gzipped)
RESULT_DIR="$PROJECT_DIR/results"               # assignments, counts, calibration
LOG_DIR="$PROJECT_DIR/logs"                     # runtime and peak-memory logs (for R1)

# Samples that are aquarium positive controls (names = barcode folder or file stem).
# Used by the calibration step. Example: AQ_SAMPLES="barcode19 barcode20 barcode21"
AQ_SAMPLES="barcode20 barcode21 barcode22"

# 2. COMPUTE -----------------------------------------------------------------
THREADS=8
NCBI_API_KEY=""        # optional; with a key, NCBI allows 10 requests/s instead of 3
                       # (free from your NCBI account settings page)

# 3. READ FILTERING ----------------------------------------------------------
# First filter: reads 7-10.5 kb with mean quality >= Q20 (BLAST_MIN_Q). Only these reads are BLASTed.
# Q20 was chosen by calibration on the aquarium samples (per-read identity by quality; see README).
BLAST_MIN_Q=20
MIN_LEN=7000
MAX_LEN=10500
FRAGMENT_SPLIT=9200    # reads < 9200 bp = 8.0 kb (12S) fragment; >= 9200 = 9.5 kb (COI)

# 4. BLAST RETRIEVAL ---------------------------------------------------------
# Deliberately loose: it only decides which hits are kept for the assignment step.
# Species-level thresholds are applied (and calibrated) in 04/05, not here.
PERC_IDENTITY=85
WORD_SIZE=64           # 9x faster than 28 with the same reads hit (200-read test);
                       # 100 (original script) risks missing seeds in raw ONT reads
EVALUE=1e-20
MAX_TARGET_SEQS=100    # 5 is too few for an LCA: one species' many identical
                       # mitogenomes can crowd out its congeners

export BLASTDB="$REF_DIR"   # lets BLAST find the taxonomy files (taxdb)
if [ -n "$NCBI_API_KEY" ]; then export NCBI_API_KEY; fi
