#!/usr/bin/env python3
"""Read accounting for the supplement: reads kept and lost at each step, per sample.

Steps (reads):
  1 basecalled_pass     MinKNOW pass reads (complete records; count_raw_reads.sh)
  2 length_q20          7-10.5 kb and mean Q >= 20 (the first filter)
  3 blast_hit           >= 1 reference with >= 2,000 bp aligned
  4 assigned            classified at family level or finer (>= 90% identity)
  5 species             classified to species (>= 98% identity, LCA within 0.5%)
  6 species_after_bleed species-level reads left after the 1% barcode-bleed filter

Usage (from the project folder, env active):
  python3 08_read_accounting.py raw_read_totals.tsv
Writes results/stats/read_accounting_by_sample.tsv and read_accounting_summary.tsv.
"""
import os, re, sys
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from lrpcr_assign import read_config  # noqa: E402

raw_path = sys.argv[1] if len(sys.argv) > 1 else "raw_read_totals.tsv"
cfg = read_config(os.path.join(HERE, "config.sh"))
res = cfg["RESULT_DIR"]
out_dir = os.path.join(res, "stats")
os.makedirs(out_dir, exist_ok=True)

bc = lambda s: (re.search(r"barcode\d+", str(s)) or [None])[0]

meta = pd.read_csv(os.path.join(HERE, "lists", "sample_metadata.csv"), encoding="utf-8-sig")
meta = meta.rename(columns={"demulti": "barcode"})[["barcode", "SampleID", "Source"]]

raw_all = pd.read_csv(raw_path, sep="\t")
raw = raw_all[raw_all["barcode"].astype(str).str.match(r"barcode\d+$")]
passed = raw[raw["folder"].str.startswith("fastq_pass")].groupby("barcode")["reads"].sum()
failed = raw[raw["folder"].str.startswith("fastq_fail")].groupby("barcode")["reads"].sum()

a = pd.read_csv(os.path.join(res, "assignment", "read_assignments.tsv.gz"), sep="\t",
                usecols=["sample", "read", "status", "rank"])
a["barcode"] = a["sample"].map(bc)
tc = pd.read_csv(os.path.join(res, "assignment", "taxon_counts.tsv"), sep="\t")
tc["barcode"] = tc["sample"].map(bc)

rows = []
for b in meta["barcode"]:
    s = a[a["barcode"] == b]
    t = tc[(tc["barcode"] == b) & (tc["rank"] == "species")]
    rows.append({
        "barcode": b,
        "basecalled_fail": int(failed.get(b, 0)),
        "basecalled_pass": int(passed.get(b, 0)),
        "length_q20": len(s),   # reads with mean Q >= 20 and 7-10.5 kb (assignment input)
        "blast_hit": int((s["status"] != "no_hit").sum()),
        "assigned": int((s["status"] == "assigned").sum()),
        "species": int((s["rank"] == "species").sum()),
        "species_after_bleed": int(t.loc[t["kept"], "reads"].sum()),
    })
df = meta.merge(pd.DataFrame(rows), on="barcode")
steps = ["basecalled_pass", "length_q20", "blast_hit", "assigned", "species", "species_after_bleed"]

def pct(d):
    d = d.copy()
    for prev, cur in zip(steps, steps[1:]):
        d[f"{cur}_pct_of_prev"] = (100 * d[cur] / d[prev]).where(d[prev] > 0).round(1)
    d["final_pct_of_pass"] = (100 * d[steps[-1]] / d[steps[0]]).where(d[steps[0]] > 0).round(1)
    return d

pct(df).to_csv(os.path.join(out_dir, "read_accounting_by_sample.tsv"), sep="\t", index=False)

grp = df.assign(group=df["Source"].where(df["Source"].isin(["Field", "Aquarium"]), "Control"))
summ = grp.groupby("group")[["basecalled_fail"] + steps].sum()
summ.loc["All samples"] = summ.sum()
summ = pct(summ.reset_index())
summ.to_csv(os.path.join(out_dir, "read_accounting_summary.tsv"), sep="\t", index=False)

# Unclassified (no barcode) reads are reported separately.
unc = raw_all[raw_all["barcode"].astype(str).str.contains("unclassified")].groupby("folder")["reads"].sum()
with pd.option_context("display.width", 200, "display.max_columns", 30):
    print(summ.T.to_string(header=False))
    if len(unc):
        print("\nUnclassified (no barcode) reads:", unc.to_dict())
print(f"\nWrote {out_dir}/read_accounting_by_sample.tsv and read_accounting_summary.tsv")
