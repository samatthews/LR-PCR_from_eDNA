#!/usr/bin/env python3
"""
Step 10 - where are tank species lost? Follows every aquarium read that could belong to each
tank species through the assignment steps, for one assay (LR-PCR or MiFish).

For each tank species S (lists/aquarium_species.tsv) and fragment, reads are counted by fate:
  no_long_alignment    the read's best hit is S, but no hit reaches the minimum aligned length
  candidate_reads      S is among the best hits (within the LCA window) of a read
    assigned_to_S        ... and the read was assigned to S (before the bleed filter)
    assigned_other_sp    ... but a local-species rule assigned it to another species
    tie_above_sp_id      ... identity >= species threshold, but other species tied -> genus/family
    below_sp_id          ... best identity below the species threshold -> genus/family
                             (split: within 1% of the threshold, and further below)
    unassigned           ... no rank could be agreed
  outscored_reads      S matches the read at >= the species threshold, but another species
                       matches better (outside the LCA window), so S is not a candidate
  bleed_removed_reads  species-level reads of S removed by the barcode-bleed filter
                       (and in how many aquarium samples)
  final_reads / final_replicates   what is reported

Also writes, for all genus-level reads in the aquarium samples, the candidate species and
fragment (genus_level_by_fragment.tsv) - e.g. which fragment the Sebastes ties come from.

Inputs: results of lrpcr_assign.py for that config (read_assignments.tsv.gz, taxon_counts.tsv,
params.json), the BLAST results and reads it used, and lists/aquarium_species.tsv.
MiFish aquarium samples = AQ_SAMPLES in config_mifish.sh; LR-PCR = Aquarium barcodes in
lists/sample_metadata.csv.

Outputs (RESULT_DIR/diagnostics/):
  tank_species_fate.tsv        one row per tank species x fragment
  tie_partners.tsv             tank species x the species it tied with, reads
  genus_level_by_fragment.tsv  genus-level reads by fragment, genus and candidate species

Usage:  python3 10_diagnose_detections.py --config config_mifish.sh
        python3 10_diagnose_detections.py --config config.sh
Time: about as long as lrpcr_assign.py (it reloads the BLAST hits from the cache).
"""
import argparse, json, os, re, sys
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import lrpcr_assign as la  # noqa: E402

KEY = ["sample", "read"]


def aquarium_samples(cfg, samples):
    if (cfg.get("ASSAY_LABEL") or "").strip():
        aq = (cfg.get("AQ_SAMPLES") or "").split()
        if not aq:
            sys.exit("Set AQ_SAMPLES in the MiFish config to the aquarium sample names.")
        return [s for s in aq if s in samples] or sys.exit("None of AQ_SAMPLES is in the reads folder.")
    meta = pd.read_csv(os.path.join(HERE, "lists", "sample_metadata.csv"), encoding="utf-8-sig")
    bc = set(meta.loc[meta["Source"].astype(str).str.strip() == "Aquarium", "demulti"].astype(str))
    return [s for s in samples if (re.search(r"barcode\d+", s) or [None])[0] in bc]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", default=os.path.join(HERE, "config.sh"))
    ap.add_argument("--tank", default=os.path.join(HERE, "lists", "aquarium_species.tsv"))
    a = ap.parse_args()

    cfg = la.read_config(a.config)
    adir = os.path.join(cfg["RESULT_DIR"], "assignment")
    p = {**la.DEFAULTS, **json.load(open(os.path.join(adir, "params.json")))}
    out = os.path.join(cfg["RESULT_DIR"], "diagnostics")
    os.makedirs(out, exist_ok=True)

    samples = aquarium_samples(cfg, la.list_samples(cfg["READ_DIR"]))
    print(f"Aquarium samples: {', '.join(samples)}", file=sys.stderr)
    print(f"Settings: min_aln {p['min_aln']}, delta {p['delta']}, species {p['sp_id']}%, "
          f"bleed {p['bleed']}", file=sys.stderr)

    tank = pd.read_csv(a.tank, sep="\t")
    tank_sp = list(tank["species"])

    # --- reads and hits, exactly as lrpcr_assign.py sees them
    lin = la.load_lineages(os.path.join(cfg["REF_DIR"], "taxid_lineage.tsv"))[["taxid", "species"]]
    hits, info = la.load_all(samples, cfg["BLAST_DIR"], cfg["READ_DIR"], os.path.join(cfg["RESULT_DIR"], "cache"))
    reads = info[(info["meanq"] >= p["min_q"]) | info["meanq"].isna()].copy()
    if p.get("assay_label"):
        reads["fragment"] = p["assay_label"]
    else:
        reads["fragment"] = np.where(reads["qlen"] < p["fragment_split"], "8.0kb", "9.5kb")
    hits = hits.merge(reads[KEY + ["fragment"]], on=KEY).merge(lin, on="taxid", how="left")

    # Reads with no hit reaching min_aln: which species was their best (longest, then best identity) hit?
    has_long = hits.loc[hits["aln"] >= p["min_aln"], KEY].drop_duplicates()
    short = hits.merge(has_long.assign(_l=1), on=KEY, how="left")
    short = short[short["_l"].isna()].sort_values(["aln", "pid"], ascending=False).drop_duplicates(KEY)
    no_long = short.groupby(["fragment", "species"]).size().rename("no_long_alignment")

    # Candidates: hits >= min_aln within delta of the read's best identity.
    h = hits[hits["aln"] >= p["min_aln"]].copy()
    h["best_pid"] = h.groupby(KEY)["pid"].transform("max")
    cand = h[h["pid"] >= h["best_pid"] - p["delta"]]
    rc = cand.dropna(subset=["species"]).groupby(KEY + ["fragment", "species"], as_index=False)["best_pid"].first()
    cset = cand.dropna(subset=["species"]).groupby(KEY)["species"].agg(lambda s: frozenset(s)).rename("cset")

    # Species that match at >= sp_id but are outscored (not candidates).
    good = h[(h["pid"] >= p["sp_id"])].dropna(subset=["species"])[KEY + ["fragment", "species"]].drop_duplicates()
    good = good.merge(rc[KEY + ["species"]].assign(_c=1), on=KEY + ["species"], how="left")
    outscored = good[good["_c"].isna()].groupby(["fragment", "species"]).size().rename("outscored_reads")

    # Fate of each candidate read, from the assignments.
    asg = pd.read_csv(os.path.join(adir, "read_assignments.tsv.gz"), sep="\t",
                      usecols=lambda c: c in {"sample", "read", "rank", "taxon", "status"}, dtype={"read": str})
    asg = asg[asg["sample"].isin(samples)]
    rc = rc.merge(asg, on=KEY, how="left").merge(cset, left_on=KEY, right_index=True, how="left")
    fate = np.select(
        [(rc["rank"] == "species") & (rc["taxon"] == rc["species"]),
         rc["rank"] == "species",
         rc["status"] != "assigned",
         rc["best_pid"] >= p["sp_id"],
         rc["best_pid"] >= p["sp_id"] - 1],
        ["assigned_to_S", "assigned_other_sp", "unassigned", "tie_above_sp_id", "below_sp_id_within_1pct"],
        default="below_sp_id_further")
    rc["fate"] = fate
    tab = rc.pivot_table(index=["fragment", "species"], columns="fate", values="read",
                         aggfunc="count", fill_value=0)
    tab["candidate_reads"] = tab.sum(axis=1)

    # Tie partners for tank species.
    ties = rc[rc["fate"] == "tie_above_sp_id"].copy()
    ties["partners"] = [";".join(sorted(cs - {sp})) for cs, sp in zip(ties["cset"], ties["species"])]
    ties = ties[ties["species"].isin(tank_sp)].groupby(["fragment", "species", "partners"]).size() \
        .rename("reads").reset_index().sort_values(["species", "fragment", "reads"], ascending=[True, True, False])
    ties.to_csv(os.path.join(out, "tie_partners.tsv"), sep="\t", index=False)

    # Bleed filter and final counts.
    tc = pd.read_csv(os.path.join(adir, "taxon_counts.tsv"), sep="\t")
    tc = tc[tc["sample"].isin(samples) & (tc["rank"] == "species")]
    tc["kept"] = tc["kept"].astype(str).str.lower().eq("true")
    g = tc.groupby(["fragment", "taxon"])
    bleed = pd.DataFrame({
        "bleed_removed_reads": g.apply(lambda d: d.loc[~d["kept"], "reads"].sum()),
        "bleed_removed_samples": g.apply(lambda d: int((~d["kept"]).sum())),
        "final_reads": g.apply(lambda d: d.loc[d["kept"], "reads"].sum()),
        "final_replicates": g.apply(lambda d: int(d["kept"].sum())),
    })
    bleed.index.names = ["fragment", "species"]

    frags = sorted(reads["fragment"].unique())
    idx = pd.MultiIndex.from_product([frags, tank_sp], names=["fragment", "species"])
    res = pd.DataFrame(index=idx).join(no_long).join(tab).join(outscored).join(bleed).fillna(0)
    cols = ["no_long_alignment", "candidate_reads", "assigned_to_S", "assigned_other_sp", "tie_above_sp_id",
            "below_sp_id_within_1pct", "below_sp_id_further", "unassigned", "outscored_reads",
            "bleed_removed_reads", "bleed_removed_samples", "final_reads", "final_replicates"]
    for c in cols:
        if c not in res:
            res[c] = 0
    res = res[cols].astype(int).reset_index()
    res = res.merge(tank[["species", "tank_census"]], on="species", how="left")
    res = res[["fragment", "species", "tank_census"] + cols].sort_values(["tank_census", "species", "fragment"],
                                                                         ascending=[False, True, True])
    res.to_csv(os.path.join(out, "tank_species_fate.tsv"), sep="\t", index=False)

    # Genus-level reads by fragment and candidate set (e.g. the Sebastes ties).
    full = pd.read_csv(os.path.join(adir, "read_assignments.tsv.gz"), sep="\t",
                       usecols=lambda c: c in {"sample", "read", "rank", "taxon", "fragment", "cand_species"},
                       dtype={"read": str})
    full = full[full["sample"].isin(samples) & (full["rank"] == "genus")]
    (full.groupby(["fragment", "taxon", "cand_species"]).size().rename("reads").reset_index()
         .sort_values(["taxon", "fragment", "reads"], ascending=[True, True, False])
         .to_csv(os.path.join(out, "genus_level_by_fragment.tsv"), sep="\t", index=False))

    # Length window (MiFish): how many raw reads fell outside MIN_LEN-MAX_LEN?
    rcp = os.path.join(cfg["RESULT_DIR"], "read_counts.tsv")
    if p.get("assay_label") and os.path.exists(rcp):
        r = pd.read_csv(rcp, sep="\t")
        r = r[r["sample"].isin(samples)]
        print(f"\nLength window: {int(r['too_short'].sum()):,} too short and {int(r['too_long'].sum()):,} too long, "
              f"of {int(r['raw_reads'].sum()):,} raw aquarium reads", file=sys.stderr)

    pd.set_option("display.width", 250)
    print(res.to_string(index=False))
    print(f"\nWrote {out}/tank_species_fate.tsv, tie_partners.tsv, genus_level_by_fragment.tsv")


if __name__ == "__main__":
    main()
