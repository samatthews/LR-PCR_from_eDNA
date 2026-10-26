# LR-PCR BLAST-only pipeline

This pipeline replaces amplicon_sorter clustering: each filtered read is classified directly by BLAST against one database built from your 12S and COI GI lists. There are no consensus sequences (so no chimeric consensus, R2 point 6) and no manual Geneious step (R1 L195–205). Thresholds are calibrated on the aquarium controls rather than chosen in advance.

## Run order (on your Mac)

| Step | Command | What it does | Time |
|---|---|---|---|
| 0 | `bash 00_setup_env.sh` then `conda activate lrpcr` | Installs BLAST, EDirect, seqkit, taxonkit, pandas | 10 min |
| – | edit `config.sh` | Set `FASTQ_DIR`, `AQ_SAMPLES`, `THREADS`, optionally `NCBI_API_KEY` | 5 min |
| 1 | `bash 01_build_refdb.sh` | Downloads ~710k records (union of both lists) straight from NCBI (Python, no EDirect), checks each chunk is complete, and builds `refdb/vert_mito` with taxonomy. Resumable. | 1–3 h, unattended |
| 2 | `bash 02_filter_reads.sh` | From the raw barcode folders: one copy of each chunk (`.fastq.gz`, else the most complete of `.fastq` / `_sanitized.fastq`), duplicate read IDs removed, then the first filter (mean Q ≥ 20, 7–10.5 kb). Writes `read_counts.tsv` and `raw_file_choice.tsv` | minutes |
| 3a | `bash 03_blast_reads.sh --aquarium --subsample 5000` | **Pilot:** BLAST 5,000 random reads per aquarium sample | ~20–35 min |
| 5 | `python3 05_calibrate.py --subsample 5000` | **Go/no-go:** grid of Q × min alignment × species identity × LCA delta, scored against the tank list | minutes |
| 6 | `python3 06_reference_coverage.py` | Mitogenome / 12S / COI records per tank species | seconds |
| 3b | `bash 03_blast_reads.sh` | BLAST all samples in full (skips finished ones) | ~6 h at word size 64; run overnight |
| 4 | `python3 lrpcr_assign.py --min-q .. --sp-id .. --delta .. --min-aln ..` | Final assignment with the chosen settings | minutes |
| 7 | `Rscript 07_stats_figures.R` | Species per sample (mean ± SE with zeros), PCR dropout (Table S1), 8.0 vs 9.5 kb overlap, controls, aquarium vs census; heatmap (Figs 6–7) and occupancy (Fig 8) figures | ~1 min |

Get the NCBI API key before step 1 if you can (NCBI account → Account settings → API Key Management). It lets the download use 4 parallel requests instead of 2.

Raw read totals for the supplement: `bash count_raw_reads.sh <run folder> > raw_read_totals.tsv` (pass and fail folders; same one-copy-per-chunk rule). Step 1 builds the reference database and does not use the reads, so it doesn't need re-running.

## Go/no-go on Sep 30

After step 5, send me `results/calibration/best_settings.txt`, `identity_by_q.tsv` and `results/reference_coverage.tsv`. **Go** if a setting detects at least the 8 tank species the submission found, with no more off-list species than before. If not, keep the clustering pipeline and use the calibration table to justify its thresholds.

## Outputs used in the revision

- `results/read_counts.tsv`: raw vs filtered reads per sample (R1 L238–240)
- `results/stats/read_accounting_*.tsv` (08_read_accounting.py): reads kept at each step, for Table S#
- `logs/compute_summary.tsv`: wall time and peak RAM per sample (R1 L81–84)
- `results/calibration/identity_by_q.tsv`: per-read identity to the true species by Q, which is the evidence for the Q and identity thresholds (R2 5; R3 L178, L348)
- `results/calibration/calibration_grid.tsv`: the sensitivity analysis (R4 L174–241)
- `results/reference_coverage.tsv`: missing references vs non-amplification (R1 major 2; R2 8)
- `results/assignment/species_by_sample.tsv` and `sample_summary.tsv`: inputs for the detection statistics and heatmap figures

## Assignment rules (lrpcr_assign.py)

1. HSPs are summed per read and reference. A read crossing the mitogenome linearization point is scored over its whole alignment, which replaces the 60% query-coverage filter (R3 L353).
2. References with fewer than `--min-aln` aligned bp (default 2,000, from calibration) are ignored. This stops a 170 bp MiFish record outranking a full mitogenome.
3. Candidates are references within `--delta` % identity of the best hit (default 0.5). The LCA is the deepest rank all candidates share.
4. The rank is capped by best identity: species ≥ `--sp-id` (98), genus ≥ `--gen-id` (95), family ≥ `--fam-id` (90); lower is unassigned.
5. The Hippoglossoides and Hippoglossus rules from the manuscript are kept.
6. Barcode-bleed filter at 1% (`--bleed 0.01`), set from the observed crosstalk (≤ 0.07% of a species' reads).

Q20 is the first filter (BLAST_MIN_Q=20 in config.sh). It was chosen by calibration: a one-off pilot BLAST of aquarium reads filtered with a looser quality setting compared per-read identity by quality (05_calibrate.py, identity_by_q.tsv). To reproduce that pilot, set BLAST_MIN_Q lower for steps 2–3a only.

## Things to check

- **Rockfish.** Twelve of the 20 tank species are *Sebastes*, and their mitogenomes differ very little. Many reads may stop at genus. The calibration reports `sebastes_genus_reads` so you can see this directly. It is likely part of the answer to R2 point 8.
- **Fragment split.** Reads < 9,200 bp are called the 8.0 kb fragment (as in the manuscript, L204–205).
- The Python steps were tested on synthetic BLAST output covering split HSPs, congeners, the flatfish rules, identity caps, no-hit and low-Q reads, short partial references and barcode bleed. The shell steps (NCBI download, BLAST) could not be run where they were written, so tell me about any error and I'll fix it.

## MiFish comparison (ONT MiFish vs ONT long-range PCR)

The MiFish reads (~170 bp 12S, ONT) go through the same database, BLAST search and assignment rules, so both assays are scored on the same reference. Settings are in `config_mifish.sh`, which reads everything else from `config.sh`. Outputs go to `mifish/`, so the LR-PCR results are not touched.

| Step | Command | What it does |
|---|---|---|
| – | edit `config_mifish.sh` | Set `MIFISH_DIR` (one file per sample) and `AQ_SAMPLES` (the aquarium MiFish samples) |
| 2m | `bash 02_prep_mifish.sh` | Keeps reads 140–220 bp (drops the ~250 bp bacterial 16S co-amplicons); FASTQ input is also filtered at Q20. Writes `mifish/results/read_counts.tsv` |
| 3m | `LRPCR_CONFIG=config_mifish.sh bash 03_blast_reads.sh` | Same BLAST search as LR-PCR, but word size 28 |
| 5m | `python3 05_calibrate.py --config config_mifish.sh` | Optional: threshold grid on the MiFish aquarium samples (min alignment 80/120/150 bp) |
| 4m | `python3 lrpcr_assign.py --config config_mifish.sh` | Assignment with minimum alignment 120 bp; all reads labelled "MiFish" |
| 9 | `python3 09_compare_assays.py` | Aquarium detections per assay (LR 8.0 kb, 9.5 kb, either, MiFish) against the tank list |

Differences from the LR-PCR settings, and why:
- **Length window 140–220 bp** instead of 7–10.5 kb (primer-trimmed MiFish-U inserts are ~160–190 bp).
- **No Q filter for FASTA input.** The MiFish files are FASTA, so there is no per-read quality; reads are used as delivered from the upstream QC. FASTQ input gets the same Q20 filter as LR-PCR.
- **Word size 28** instead of 64; a 64 bp exact seed is over a third of a MiFish read.
- **Minimum alignment 120 bp** instead of 2,000; 2,000 would reject every MiFish hit.
- Species identity (98%), LCA delta (0.5%), the local-species rules and the 1% bleed filter are unchanged. On a 170 bp read one mismatch is 0.6%, so a 0.5% delta keeps only references with the fewest mismatches.

`collab_blast.sh` runs step 3 on another computer: send it the step-2 FASTAs (`reads/` or `mifish/reads/`) and the `refdb/` folder; the `.blast.tsv.gz` files it returns go into `blast/` or `mifish/blast/`. Use `WORD_SIZE=28` for MiFish.

## BLAST on another computer or database (lab server)

1. `bash 02_filter_reads.sh` here, then copy `reads/*.fasta` to the server.
2. On the server: `DB=/path/to/db_prefix THREADS=32 bash collab_blast.sh fasta_folder/ blast_out/` (MiFish: add `WORD_SIZE=28`). Same search settings as step 3.
3. Copy `blast_out/*.blast.tsv.gz` back into an empty `blast/` (move any earlier vert_mito results to another folder first, so the two databases are never mixed). Keep `blast_version.txt` and `compute_summary.tsv` for the Methods.
4. If the database is not `refdb/vert_mito`: `bash 04a_lineages_from_blast.sh`, which adds lineages for the new taxids so the assignment step can classify them.
5. Re-run the calibration (`05_calibrate.py`) before step 4: the thresholds were chosen with vert_mito and may shift with a larger database. `06_reference_coverage.py` still describes vert_mito, not the server database.
