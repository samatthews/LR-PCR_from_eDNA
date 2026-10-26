#!/usr/bin/env python3
"""
Step 11 - barcode-bleed sensitivity table (supplementary table; reviewer concern that
a fixed cutoff is arbitrary).

taxon_counts.tsv from lrpcr_assign.py keeps every sample x fragment x taxon row,
including those removed by the bleed filter, with each row's share of that taxon's
reads (share_of_taxon). This script re-applies the filter at a range of cutoffs to
that file, so nothing needs re-running: BLAST and the assignment settings stay as they
are, and only the bleed cutoff varies. It counts species the same way as
07_stats_figures.R (species-level calls, either fragment, zeros included).

Inputs:
  RESULT_DIR/assignment/taxon_counts.tsv   (from lrpcr_assign.py, current settings)
  lists/sample_metadata.csv                (Field / Aquarium / NC)
  lists/aquarium_species.tsv               (tank census)

Output: RESULT_DIR/stats/bleed_sensitivity.tsv, one row per cutoff:
  field species per sample (mean, SE, range), total field species, aquarium species per
  sample, total aquarium species, tank species detected (census > 0), off-list species in
  the aquarium, and detections supported by 1 or 2 reads (the low-level cross-sample reads
  the filter is meant to remove).

Usage:  python3 11_bleed_sensitivity.py
        python3 11_bleed_sensitivity.py --cutoffs 0 0.001 0.005 0.01 0.03
"""
import argparse, os, sys
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import lrpcr_assign as la  # noqa: E402

DEFAULT_CUTOFFS = [0, 0.001, 0.0025, 0.005, 0.01, 0.02, 0.03]


def se(x):
    x = np.asarray(x, dtype=float)
    return x.std(ddof=1) / np.sqrt(len(x)) if len(x) > 1 else np.nan


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", default=os.path.join(HERE, "config.sh"))
    ap.add_argument("--cutoffs", nargs="+", type=float, default=DEFAULT_CUTOFFS,
                    help="bleed cutoffs as fractions (0.005 = 0.5%%)")
    a = ap.parse_args()

    cfg = la.read_config(a.config)
    path = os.path.join(cfg["RESULT_DIR"], "assignment", "taxon_counts.tsv")
    if not os.path.exists(path):
        sys.exit(f"No {path} - run lrpcr_assign.py first")
    c = pd.read_csv(path, sep="\t")
    c = c[c["rank"] == "species"].copy()
    c["barcode"] = c["sample"].astype(str).str.extract(r"(barcode\d+)", expand=False)

    meta = pd.read_csv(os.path.join(HERE, "lists", "sample_metadata.csv"), encoding="utf-8-sig")
    meta = meta.rename(columns={"demulti": "barcode"})[["barcode", "SampleID", "Source"]]
    c = c.merge(meta, on="barcode", how="inner")
    tank = pd.read_csv(os.path.join(HERE, "lists", "aquarium_species.tsv"), sep="\t")
    in_tank = set(tank.loc[tank["tank_census"] > 0, "species"])
    on_list = set(tank["species"])

    field_bc = meta.loc[meta["Source"] == "Field", "barcode"].tolist()
    aq_bc = meta.loc[meta["Source"] == "Aquarium", "barcode"].tolist()

    rows = []
    for cut in sorted(a.cutoffs):
        k = c[c["share_of_taxon"] >= cut]
        # either fragment: a species counts once per sample
        det = k.groupby(["barcode", "Source", "taxon"], as_index=False)["reads"].sum()
        per = det.groupby("barcode")["taxon"].nunique()
        f = per.reindex(field_bc, fill_value=0)
        q = per.reindex(aq_bc, fill_value=0)
        fd, ad = det[det["Source"] == "Field"], det[det["Source"] == "Aquarium"]
        aq_sp = set(ad["taxon"])
        rows.append({
            "bleed_cutoff_pct": round(cut * 100, 3),
            "field_species_per_sample_mean": round(f.mean(), 2),
            "field_species_per_sample_se": round(se(f), 2),
            "field_species_per_sample_range": f"{f.min()}-{f.max()}",
            "field_species_total": fd["taxon"].nunique(),
            "field_detections": len(fd),
            "field_detections_1_or_2_reads": int((fd["reads"] <= 2).sum()),
            "aquarium_species_per_sample_mean": round(q.mean(), 2),
            "aquarium_species_per_sample_se": round(se(q), 2),
            "aquarium_species_total": len(aq_sp),
            "tank_species_detected": f"{len(aq_sp & in_tank)} of {len(in_tank)}",
            "aquarium_offlist_species": len(aq_sp - on_list),
            "aquarium_detections_1_or_2_reads": int((ad["reads"] <= 2).sum()),
        })
    out = pd.DataFrame(rows)
    od = os.path.join(cfg["RESULT_DIR"], "stats")
    os.makedirs(od, exist_ok=True)
    out.to_csv(os.path.join(od, "bleed_sensitivity.tsv"), sep="\t", index=False)
    print(out.to_string(index=False))
    print(f"\nCurrent default cutoff: {la.DEFAULTS['bleed'] * 100:g}%")
    print(f"Wrote {os.path.join(od, 'bleed_sensitivity.tsv')}")


if __name__ == "__main__":
    main()
