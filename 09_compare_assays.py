#!/usr/bin/env python3
"""
Step 9 - compare the aquarium (WOWW tank) detections of the two ONT assays scored
against the same reference database and assignment rules:
  LR-PCR 8.0 kb, LR-PCR 9.5 kb, LR-PCR either fragment, and MiFish (~170 bp 12S).

Inputs (run lrpcr_assign.py for both configs first):
  results/assignment/taxon_counts.tsv          (LR-PCR, config.sh)
  mifish/results/assignment/taxon_counts.tsv   (MiFish, config_mifish.sh)
  lists/aquarium_species.tsv                   (tank list and census)
LR-PCR aquarium samples are the Aquarium barcodes in lists/sample_metadata.csv;
MiFish aquarium samples are AQ_SAMPLES in config_mifish.sh.

Only species-level detections that pass the barcode-bleed filter are counted, as in
lrpcr_assign.py.

Outputs (mifish/results/comparison/):
  aquarium_species_by_assay.tsv  every tank species plus any off-list species: replicates
                                 detected and reads, per assay
  aquarium_assay_summary.tsv     per assay: replicates, reads per replicate, tank species
                                 detected (total and mean +/- SE per replicate), off-list species

Usage:  python3 09_compare_assays.py
        python3 09_compare_assays.py --lr-config config.sh --mifish-config config_mifish.sh
"""
import argparse, os, re, sys
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import lrpcr_assign as la  # noqa: E402


def species_counts(cfg_path, label):
    cfg = la.read_config(cfg_path)
    path = os.path.join(cfg["RESULT_DIR"], "assignment", "taxon_counts.tsv")
    if not os.path.exists(path):
        sys.exit(f"No {path} - run: python3 lrpcr_assign.py --config {os.path.basename(cfg_path)}")
    c = pd.read_csv(path, sep="\t")
    c = c[c["kept"].astype(str).str.lower().eq("true") & (c["rank"] == "species")]
    return cfg, c


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--lr-config", default=os.path.join(HERE, "config.sh"))
    ap.add_argument("--mifish-config", default=os.path.join(HERE, "config_mifish.sh"))
    ap.add_argument("--tank", default=os.path.join(HERE, "lists", "aquarium_species.tsv"))
    a = ap.parse_args()

    lr_cfg, lr = species_counts(a.lr_config, "LR-PCR")
    mf_cfg, mf = species_counts(a.mifish_config, "MiFish")

    # Aquarium samples for each assay.
    meta = pd.read_csv(os.path.join(HERE, "lists", "sample_metadata.csv"), encoding="utf-8-sig")
    aq_bc = set(meta.loc[meta["Source"].astype(str).str.strip() == "Aquarium", "demulti"].astype(str))
    bc = lambda s: (re.search(r"barcode\d+", str(s)) or [None])[0]
    lr_all = sorted(set(lr["sample"]))
    lr_aq = sorted(s for s in lr_all if bc(s) in aq_bc)
    mf_aq = (mf_cfg.get("AQ_SAMPLES") or "").split()
    if not lr_aq:
        sys.exit("No LR-PCR aquarium samples found in taxon_counts.tsv (expected barcode20-22).")
    if not mf_aq:
        sys.exit("Set AQ_SAMPLES in config_mifish.sh to the MiFish aquarium sample names.")
    missing = [s for s in mf_aq if s not in set(mf["sample"])]
    if missing:
        print(f"Note: no species-level detections for MiFish sample(s): {', '.join(missing)}", file=sys.stderr)

    lr = lr[lr["sample"].isin(lr_aq)]
    mf = mf[mf["sample"].isin(mf_aq)]
    either = (lr.groupby(["sample", "taxon"], as_index=False)["reads"].sum().assign(fragment="LR either"))
    assays = {
        "LR 8.0 kb": (lr[lr["fragment"] == "8.0kb"], lr_aq),
        "LR 9.5 kb": (lr[lr["fragment"] == "9.5kb"], lr_aq),
        "LR either": (either, lr_aq),
        "MiFish": (mf, mf_aq),
    }

    tank = pd.read_csv(a.tank, sep="\t")
    on_tank = set(tank.loc[tank["tank_census"] > 0, "species"])
    zero = set(tank.loc[tank["tank_census"] == 0, "species"])
    detected_any = set().union(*(set(d["taxon"]) for d, _ in assays.values()))
    off_list = sorted(detected_any - on_tank - zero)

    # Species x assay table.
    rows = tank[["species", "common_name", "tank_census"]].copy()
    rows["on_list"] = "Y"
    if off_list:
        rows = pd.concat([rows, pd.DataFrame({"species": off_list, "common_name": "", "tank_census": np.nan,
                                               "on_list": "N"})], ignore_index=True)
    for name, (d, samples) in assays.items():
        reps = d.groupby("taxon")["sample"].nunique()
        reads = d.groupby("taxon")["reads"].sum()
        rows[f"{name}: replicates (of {len(samples)})"] = rows["species"].map(reps).fillna(0).astype(int)
        rows[f"{name}: reads"] = rows["species"].map(reads).fillna(0).astype(int)

    # Summary per assay.
    summ = []
    for name, (d, samples) in assays.items():
        per_rep = [len(set(d.loc[d["sample"] == s, "taxon"]) & on_tank) for s in samples]
        reads_rep = [int(d.loc[d["sample"] == s, "reads"].sum()) for s in samples]
        det = set(d["taxon"])
        n = len(samples)
        summ.append({
            "assay": name,
            "replicates": n,
            "species_level_reads_per_replicate": f"{np.mean(reads_rep):,.0f} ({min(reads_rep):,}-{max(reads_rep):,})",
            "tank_species_detected_total": f"{len(det & on_tank)} of {len(on_tank)}",
            "tank_species_per_replicate_mean_se": f"{np.mean(per_rep):.1f} ± {np.std(per_rep, ddof=1) / np.sqrt(n):.1f}"
            if n > 1 else f"{per_rep[0]}",
            "listed_census0_detected": len(det & zero),
            "off_list_detected": len(det - on_tank - zero),
            "off_list_species": "; ".join(sorted(det - on_tank - zero)),
            "tank_species_missed": "; ".join(sorted(on_tank - det)),
        })
    summ = pd.DataFrame(summ)

    out = os.path.join(mf_cfg["RESULT_DIR"], "comparison")
    os.makedirs(out, exist_ok=True)
    rows.sort_values(["on_list", "tank_census", "species"], ascending=[False, False, True]) \
        .to_csv(os.path.join(out, "aquarium_species_by_assay.tsv"), sep="\t", index=False)
    summ.to_csv(os.path.join(out, "aquarium_assay_summary.tsv"), sep="\t", index=False)
    print(summ.drop(columns=["off_list_species", "tank_species_missed"]).to_string(index=False))
    print(f"\nLR-PCR aquarium samples: {', '.join(lr_aq)}\nMiFish aquarium samples: {', '.join(mf_aq)}")
    print(f"Wrote {out}")
    print("Note: read depth differs between the assays; compare species counts with the reads per replicate in mind.")


if __name__ == "__main__":
    main()
