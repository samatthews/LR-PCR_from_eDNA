# lrpcr_blast

Code for:

> Stephanie A. Matthews, Olivia M. Scott, Yin Cheong Aden Ip, Elizabeth A. Allan, Ryan P. Kelly. Long-range PCR amplification and nanopore sequencing of 8–10 kb mitochondrial fragments from environmental DNA. *Molecular Ecology Resources* (in review).

Environmental DNA from Puget Sound field samples and an aquarium tank was amplified as two half-mitogenome fragments: 8.0 kb (with 12S) and 9.5 kb (with COI). The amplicons were sequenced on Oxford Nanopore. This pipeline classifies every read directly with BLAST against a vertebrate mitochondrial database.The scripts then produce the detection statistics, figures and supplementary tables in the paper.

The same pipeline also scores ONT MiFish (12S, ~170 bp) reads from the same aquarium tank (Ip et al. 2026), which gives a comparison between the two assays.

## Requirements

- **Command line:** BLAST+ (2.12.0 for the published run; ≥ 2.12 is fine), seqkit and taxonkit. `00_setup_env.sh` installs them into a conda environment called `lrpcr`.
- **Python:** Python 3.11 with pandas and numpy.
- **R:** R ≥ 4.3 with the tidyverse (dplyr, tidyr, readr, stringr, ggplot2) and patchwork. Fig. 1 also needs sf, marmap and ggrepel.
- **Cluster (optional):** the `*.slurm` scripts run steps 3 and 6 on a SLURM cluster. Their account, partition, e-mail and paths are set for the University of Washington Hyak cluster, so edit them before use.

```bash
bash 00_setup_env.sh
conda activate lrpcr
```

## Data

- **Reads:** the raw nanopore reads (pass and fail, all barcodes) are in the NCBI Sequence Read Archive under BioProject PRJNA1542106.
- **Database:** reads were searched against `vert_mtDNA`, which holds all vertebrate mitochondrial records in NCBI core_nt (downloaded 7 October 2026; 2,239,547 records). The  accession of every record are listed in Appendix 4 of the paper (`run_appendix4_gi_list.slurm`).
- **Sample tables:** these are in `lists/`.

| File | Contents |
|---|---|
| `sample_metadata.csv` | Barcode → sample, station, source (Field / Aquarium / negative control), basin, depth, coordinates |
| `pcr_success_tableS1.csv` | All 72 field samples: volume filtered, DNA concentration, gel band per fragment (Table S1) |
| `aquarium_species.tsv` | Tank census at sampling (Table S2) |
| `regional_species.txt` | Fishes recorded in the Salish Sea (Pietsch & Orr 2019); species-level calls outside this list are reported at genus |
| `WOAC_cruise_plan.csv` | Station coordinates for Fig. 1 *(add before running step 12)* |
| `vert_12s_target_gilist.txt`, `vert_COI_target_gilist.txt` | GI lists for the optional smaller local database (`01_build_refdb.sh`); not used for the published results |

## Configuration

All settings are in `config.sh`: paths, aquarium sample names, read filters, BLAST settings and assignment thresholds. MiFish settings are in `config_mifish.sh`, which reads everything else from `config.sh`. To use a different config, pass `--config` to the Python steps, or set `LRPCR_CONFIG=...` for the shell steps.

## Running the pipeline

| Step | Command | What it does |
|---|---|---|
| 1 | `bash 01_build_refdb.sh` | *Optional.* Builds a smaller local database (`refdb/vert_mito`) from the GI lists in `lists/`. Skip this step if you have `vert_mtDNA` or core_nt. |
| 2 | `bash 02_filter_reads.sh` | Reads each MinKNOW barcode folder, uses one copy of each chunk, removes duplicate read IDs, and keeps reads of 7.0–10.5 kb with mean Q ≥ 20. Writes `reads/<sample>.fasta` and per-read length and Q. |
| 3 | `bash 03_blast_reads.sh` | BLASTs every filtered read (word size 64, ≥ 85% identity, 100 targets). Logs runtime and peak memory. On a cluster, use `sbatch --array=1-N slurm_blast_array.sh` (which runs `collab_blast.sh`) instead. |
| 4a | `bash 04a_lineages_from_blast.sh` | Adds NCBI lineages for every taxid in the BLAST output. Needed when the database is not `refdb/vert_mito`. |
| 5 | `python3 05_calibrate.py` | Runs the assignment over a grid of minimum aligned length × species identity × LCA delta. Each setting is scored against the aquarium tank census. Also writes per-read identity by quality score. |
| 4 | `python3 lrpcr_assign.py` | Final per-read assignment and per-sample read counts, with the settings below. |
| 6 | `python3 06_reference_coverage.py --metadata ref_metadata.tsv` | Counts reference records per tank species: complete mitogenomes and other records ≥ 2 kb. Cluster version: `run_06_reference_coverage.slurm`. |
| 7 | `Rscript 07_stats_figures.R` | Species per sample, PCR dropout, overlap between the 8.0 and 9.5 kb fragments, controls, aquarium detections vs census, and the heatmaps (Figs 6–7). |
| 8 | `python3 08_read_accounting.py` | Reads kept and lost at each step, per sample (Table S3). Raw totals come from `count_raw_reads.sh`. |
| 9 | `python3 09_compare_assays.py` | Aquarium detections by LR-PCR and MiFish against the tank list (Table S7). |
| 10 | `python3 10_diagnose_detections.py` | Diagnostics: where each tank species' reads end up (identity, ties, genus-level calls). |
| 11 | `Rscript 11_fig5_species_per_sample.R` | Fig. 5: species per replicate, by fragment and location. |
| 12 | `Rscript 12_fig1_sampling_map.R` | Fig. 1: station map with West Coast inset. Downloads NOAA bathymetry on the first run. |
| – | `sbatch run_appendix4_gi_list.slurm` | Appendix 4: GI, accession, taxid and length of every database record. |

Run steps 2–5 first and step 4 after calibration. The step numbers follow the order the scripts were written. For a run on another computer or cluster, copy `reads/*.fasta` from step 2, run `collab_blast.sh` (or `slurm_blast_array.sh`) there, and put the returned `*.blast.tsv.gz` files into `blast/` before step 4a. Use one database for every sample and both assays.

### MiFish reads

```bash
bash 02_prep_mifish.sh                                   # keep 140–220 bp reads
LRPCR_CONFIG=config_mifish.sh bash 03_blast_reads.sh     # word size 28
python3 lrpcr_assign.py --config config_mifish.sh        # min. alignment 120 bp, species ≥ 99%
python3 09_compare_assays.py
```

## Assignment rules (`lrpcr_assign.py`)

1. **Combine alignments.** BLAST alignment segments (HSPs) are combined for each read–reference pair. A read that crosses the point where a circular mitogenome was linearized is therefore scored over its full alignment. No query-coverage filter is used.
2. **Minimum alignment.** References with less than 2,000 aligned bp are ignored (MiFish: 120 bp). This stops short partial records from outranking full mitogenomes.
3. **Candidates and LCA.** Candidates are the references within 0.5% identity of the best hit. The read is assigned to the lowest common ancestor of those candidates.
4. **Rank caps.** The best identity limits the rank: species ≥ 98% (MiFish 99%), genus ≥ 95%, family ≥ 90%. Reads below 90% are unassigned.
5. **Local rules.** Rules for *Hippoglossoides* and *Hippoglossus* are applied. Sand lances (*Ammodytes*, *Hyperoplus*) are reported at family (Ammodytidae).
6. **Regional list.** Species-level calls for species not on `lists/regional_species.txt` are reported at genus. Turn this off with `--no-regional`.
7. **Fragment.** Reads shorter than 9,200 bp are counted as the 8.0 kb fragment; the rest are the 9.5 kb fragment.
8. **Barcode bleed.** A taxon with less than 0.5% of its total reads in a sample is removed from that sample (`--bleed 0.005`).

Every threshold is a command-line option (`python3 lrpcr_assign.py --help`). The calibration grid in step 5 shows how the results change across settings.

## Main outputs

| Output | Paper |
|---|---|
| `results/assignment/read_assignments.tsv.gz`, `taxon_counts.tsv`, `species_by_sample.tsv` | Per-read and per-sample detections |
| `results/stats/read_accounting_*.tsv` | Table S3 |
| `results/calibration/identity_by_q.tsv` | Table S4 |
| `results/calibration/calibration_grid.tsv` | Table S5 |
| `reference_coverage.tsv` | Table S6 |
| `mifish/results/comparison/` | Table S7 |
| `results/stats/*.csv`, `stats_summary.txt` | Results text (species per sample, PCR dropout, fragment overlap) |
| `results/figures/fig1_sampling_map.*` | Fig. 1 |
| `results/figures/fig5_species_per_sample.*` | Fig. 5 |
| `results/figures/fig6_heatmap_field.*`, `fig7_heatmap_aquarium.*` | Figs 6–7 |
| `logs/compute_summary.tsv` | Runtime and memory (Methods) |

## Citation

If you use this code, please cite the paper above. 
