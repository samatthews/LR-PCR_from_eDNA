#!/bin/bash
# ---------------------------------------------------------------------------
# MiFish (ONT, ~170 bp 12S) settings. Everything not set here comes from
# config.sh, so both assays use the same reference database (refdb/vert_mito),
# taxonomy and assignment rules. Outputs go to mifish/ so nothing from the
# LR-PCR run is overwritten.
#
# Run order (from this folder, conda env active):
#   bash 02_prep_mifish.sh
#   LRPCR_CONFIG=config_mifish.sh bash 03_blast_reads.sh
#   python3 05_calibrate.py --config config_mifish.sh      (optional threshold check)
#   python3 lrpcr_assign.py --config config_mifish.sh
#   python3 09_compare_assays.py
# ---------------------------------------------------------------------------
source "$(dirname "${BASH_SOURCE[0]}")/config.sh"

ASSAY_LABEL="MiFish"        # one amplicon: all reads get this label instead of 8.0kb/9.5kb

# 1. PATHS -------------------------------------------------------------------
MIFISH_DIR="/Users/stephaniematthews/Documents/Rockfish/MinION_steph/MolEcolRevisions/lrpcr_blast/mifish/demux/ontbarcoder"   # one file per sample: .fasta/.fa(.gz), or .fastq(.gz)
                                     # the file name (minus extension) becomes the sample name
MIFISH_ROOT="$PROJECT_DIR/mifish"
READ_DIR="$MIFISH_ROOT/reads"
BLAST_DIR="$MIFISH_ROOT/blast"
RESULT_DIR="$MIFISH_ROOT/results"
LOG_DIR="$MIFISH_ROOT/logs"
# REF_DIR is inherited from config.sh: the same database as the LR-PCR reads.

# Aquarium (WOWW tank) MiFish samples, by sample name (file name without extension).
# Example: AQ_SAMPLES="18_mifish_bt1_a_all 19_mifish_bt1_b_all 20_mifish_bt1_c_all"
AQ_SAMPLES="18_mifish_bt1_a_all 18_mifish_bt1_b_all 18_mifish_bt1_c_all 19_mifish_bt2_a_all 19_mifish_bt2_b_all 19_mifish_bt2_c_all 20_mifish_bt3_a_all 20_mifish_bt3_b_all 20_mifish_bt3_c_all"

# 2. READ FILTERING ----------------------------------------------------------
# Primer-trimmed MiFish-U inserts are ~160-190 bp. The window drops the ~250 bp
# bacterial 16S co-amplicons (about a third of reads in the example file).
MIN_LEN=140
MAX_LEN=220
# BLAST_MIN_Q (Q20, from config.sh) is applied only to FASTQ input. FASTA input has no
# per-read quality, so its reads are kept as delivered (QC done in the upstream pipeline).
FRAGMENT_SPLIT=0            # not used for MiFish

# 3. BLAST RETRIEVAL ---------------------------------------------------------
# Same as LR-PCR except word size: a 64 bp exact seed is over a third of a 170 bp read.
WORD_SIZE=28                # megablast default

# 4. ASSIGNMENT --------------------------------------------------------------
# Minimum aligned bp per reference. The LR-PCR value (2,000) would reject every MiFish hit.
ASSIGN_MIN_ALN=120
