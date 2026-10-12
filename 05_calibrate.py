#!/usr/bin/env python3
"""
Step 5 - calibrate thresholds on the aquarium positive controls (the go/no-go test).

Runs the assignment over a grid of settings for the AQ_SAMPLES in config.sh and scores
each against the tank inventory (lists/aquarium_species.tsv):
  tank_detected       species with tank census > 0 detected (of 20)
  listed_zero_census  species on the Aquarium list with census 0 detected (of 5)
  off_list            species detected that are not on the list at all (false-positive proxy)
  off_list_read_share share of species-level reads going to off-list species
  species_read_share  share of reads resolved to species
  sebastes_genus_reads reads stopping at genus Sebastes (rockfish mitogenomes are very similar)

Also writes the per-read identity distribution to tank species by mean-Q bin, which is
the evidence for where the species threshold should sit (R2 5, R3 L178/L348).

Usage:  python3 05_calibrate.py --config config_mifish.sh   (MiFish aquarium samples)
        python3 05_calibrate.py --subsample 5000   (after 03_blast_reads.sh --aquarium --subsample 5000)
        python3 05_calibrate.py            (uses AQ_SAMPLES from config.sh, full samples)
        python3 05_calibrate.py --quick     (smaller grid, ~4x faster)
Outputs: results/calibration/calibration_grid.tsv, identity_by_q.tsv, best_settings.txt
"""
import argparse, itertools, os, sys
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lrpcr_assign as la

HERE = os.path.dirname(os.path.abspath(__file__))


def score(assigned, counts, tank):
    on_tank = set(tank.loc[tank["tank_census"] > 0, "species"])
    zero = set(tank.loc[tank["tank_census"] == 0, "species"])
    listed = on_tank | zero
    sp = counts[counts["kept"] & (counts["rank"] == "species")]
    det = set(sp["taxon"])
    sp_reads = sp["reads"].sum()
    off = det - listed
    n = len(assigned)
    return dict(
        reads=n,
        assigned_share=round((assigned["status"] == "assigned").sum() / n, 3) if n else np.nan,
        species_read_share=round((assigned["rank"] == "species").sum() / n, 3) if n else np.nan,
        tank_detected=len(det & on_tank),
        listed_zero_census=len(det & zero),
        off_list=len(off),
        off_list_read_share=round(sp.loc[sp["taxon"].isin(off), "reads"].sum() / sp_reads, 3) if sp_reads else np.nan,
        sebastes_genus_reads=int(((assigned["rank"] == "genus") & (assigned["taxon"] == "Sebastes")).sum()),
        off_list_species="; ".join(sorted(off)),
        tank_species_missed="; ".join(sorted(on_tank - det)),
    )


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", default=os.path.join(HERE, "config.sh"))
    ap.add_argument("--samples", nargs="*", help="default: AQ_SAMPLES from config.sh")
    ap.add_argument("--tank", default=os.path.join(HERE, "lists", "aquarium_species.tsv"))
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--subsample", type=int, help="use the reads BLASTed with 03_blast_reads.sh --subsample N")
    a = ap.parse_args()

    cfg = la.read_config(a.config)
    samples = a.samples or cfg.get("AQ_SAMPLES", "").split()
    if not samples:
        sys.exit("Set AQ_SAMPLES in config.sh (or pass --samples) to the aquarium control samples.")
    read_dir, blast_dir, cache = cfg["READ_DIR"], cfg["BLAST_DIR"], os.path.join(cfg["RESULT_DIR"], "cache")
    out = os.path.join(cfg["RESULT_DIR"], "calibration")
    if a.subsample:
        sub = os.path.join(os.path.dirname(cfg["READ_DIR"].rstrip("/")), f"subsample_{a.subsample}")
        read_dir, blast_dir = os.path.join(sub, "reads"), os.path.join(sub, "blast")
        cache, out = os.path.join(sub, "cache"), os.path.join(sub, "calibration")
        if not os.path.isdir(blast_dir):
            sys.exit(f"No {blast_dir} - run 03_blast_reads.sh --aquarium --subsample {a.subsample} "
                     f"(with LRPCR_CONFIG set to this config) first")
    os.makedirs(out, exist_ok=True)
    tank = pd.read_csv(a.tank, sep="\t")
    lineages = la.load_lineages(os.path.join(cfg["REF_DIR"], "taxid_lineage.tsv"))
    try:
        hits, info = la.load_all(samples, blast_dir, read_dir, cache)
    except FileNotFoundError as e:
        sys.exit(str(e))
    split = int(cfg.get("FRAGMENT_SPLIT") or 9200)
    assay = (cfg.get("ASSAY_LABEL") or "").strip()
    base_aln = la.min_aln_from(cfg)
    no_q = info["meanq"].isna().all()
    print(f"{len(info)} reads, {hits['read'].nunique()} with hits, in {len(samples)} aquarium samples")

    # 1. Identity of each read's best hit, for reads whose best hit is a tank species.
    lin = lineages[["taxid", "species"]]
    h = hits[hits["aln"] >= base_aln]
    best = h.loc[h.groupby(["sample", "read"])["pid"].idxmax()].merge(lin, on="taxid", how="left")
    best = best.merge(info, on=["sample", "read"])
    best = best[best["species"].isin(tank.loc[tank["tank_census"] > 0, "species"])]
    if no_q:   # FASTA input (MiFish): no per-read quality, so one row for all reads
        ident = best["pid"].describe(percentiles=[.05, .1, .25, .5]).round(2).to_frame("all reads (no Q)").T
    else:
        best["q_bin"] = pd.cut(best["meanq"], [0, 15, 18, 20, 25, 99], right=False,
                               labels=["<15", "15-18", "18-20", "20-25", ">=25"])
        ident = best.groupby("q_bin", observed=True)["pid"].describe(percentiles=[.05, .1, .25, .5]).round(2)
    ident.to_csv(os.path.join(out, "identity_by_q.tsv"), sep="\t")
    print("\nBest-hit identity (%) of reads matching tank species, by read mean Q:")
    print(ident.to_string())

    # 2. Grid.
    if assay:   # ~170 bp MiFish reads: alignments can't reach the LR-PCR minimums
        aln_grid = [120] if a.quick else [80, 120, 150]
    else:
        aln_grid = [500, 2000] if a.quick else [150, 500, 2000]
    grid = dict(min_q=[0] if no_q else [15, 18, 20], min_aln=aln_grid,
                sp_id=[96, 97, 98, 99] if not a.quick else [97, 98], delta=[0.5] if a.quick else [0.25, 0.5, 1.0])
    rows = []
    combos = list(itertools.product(*grid.values()))
    for i, vals in enumerate(combos, 1):
        p = dict(zip(grid, vals)); p["fragment_split"] = split; p["assay_label"] = assay
        p["gen_id"] = min(la.DEFAULTS["gen_id"], p["sp_id"] - 1)
        asg = la.assign(hits, info, lineages, p)
        cnt = la.count_taxa(asg)
        for frag in ["all"] + ([assay] if assay else ["8.0kb", "9.5kb"]):
            aa = asg if frag == "all" else asg[asg["fragment"] == frag]
            cc = cnt if frag == "all" else cnt[cnt["fragment"] == frag]
            rows.append({**p, "fragment": frag, **score(aa, cc, tank)})
        print(f"\r  grid {i}/{len(combos)}", end="", flush=True)
    print()
    res = pd.DataFrame(rows)
    res.to_csv(os.path.join(out, "calibration_grid.tsv"), sep="\t", index=False)

    # 3. Suggest: most tank species with no off-list species; ties -> more reads kept.
    allf = res[res["fragment"] == "all"].sort_values(
        ["off_list", "tank_detected", "off_list_read_share", "reads"], ascending=[True, False, True, False])
    with open(os.path.join(out, "best_settings.txt"), "w") as f:
        f.write("Top settings (fewest off-list species, then most tank species, then most reads):\n")
        f.write(allf.head(10).drop(columns=["off_list_species", "tank_species_missed", "fragment_split",
                                            "assay_label"]).to_string(index=False) + "\n")
    print(open(os.path.join(out, "best_settings.txt")).read())
    if not assay:
        print("Submitted pipeline for comparison: 8 of 20 tank species, 5 listed-but-census-0 species.")


if __name__ == "__main__":
    main()
