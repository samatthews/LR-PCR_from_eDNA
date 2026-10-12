#!/usr/bin/env python3
"""
Step 4 - assign taxonomy to each read from its BLAST hits, then count reads per
species per sample. Also imported by 05_calibrate.py and 06_reference_coverage.py.

Per read:
  1. HSPs are summed per reference (aligned length, identity weighted by length), so a
     read that spans the artificial linearization point of a mitogenome (two HSPs)
     is scored over its full alignment rather than needing a 60% query-coverage cut.
  2. References with less than --min-aln aligned bp are ignored.
  3. Candidates = references within --delta % identity of the best hit.
  4. LCA = deepest rank on which all candidates agree.
  5. Rank is capped by the best identity: species needs >= --sp-id, genus >= --gen-id,
     family >= --fam-id; below --fam-id the read is 'unassigned'.
  6. Hippoglossoides / Hippoglossus rules from the original manuscript are applied.
Per sample: reads counted per fragment (length split) and taxon; taxa with < --bleed
of their total reads in a sample are removed from that sample (barcode bleed).

Usage (after steps 1-3):
  conda activate lrpcr
  python3 lrpcr_assign.py --config config.sh                 # defaults
  python3 lrpcr_assign.py --config config.sh --min-q 20 --sp-id 98.5 --out results/q20
  python3 lrpcr_assign.py --config config_mifish.sh          # MiFish reads (min-aln 120, one fragment)
"""
import argparse, glob, gzip, json, os, re, shlex, subprocess, sys
import numpy as np
import pandas as pd

RANKS = ["kingdom", "phylum", "class", "order", "family", "genus", "species"]

# Geographic / taxonomic rules. Where the only competing species do not occur in the
# NE Pacific (or are synonyms), assign the local species. The first two are from the
# submitted manuscript (L226-232); the rest were added in revision (Sep 30) after the
# genus-level candidate check (genus_level_candidates.tsv).
LOCAL = {
    "Hippoglossoides robustus": "Hippoglossoides elassodon",   # junior synonym (manuscript)
    "Hippoglossoides dubius": "Hippoglossoides elassodon",     # NW Pacific relative
    "Hippoglossus hippoglossus": "Hippoglossus stenolepis",    # Atlantic halibut (manuscript)
    "Merluccius angustimanus": "Merluccius productus",         # Panama hake, tropical E Pacific
    "Gadus finnmarchicus": "Gadus chalcogrammus",              # Norwegian pollock, Atlantic (and see G. morhua below)
    "Sebastes eos": "Sebastes ruberrimus",                     # pink rockfish, southern range
    "Clupea harengus": "Clupea pallasii",                      # Atlantic herring
}
# Species-level calls to a non-local species are remapped (except Atlantic halibut, which the
# manuscript only merged at genus level).
SPECIES_REMAP = {k: v for k, v in LOCAL.items() if k != "Hippoglossus hippoglossus"}
# Genus-level reads whose candidate species are only a local species plus its non-local
# relatives are assigned to the local species.
GENUS_RULES = {}
for local in set(LOCAL.values()):
    GENUS_RULES[frozenset({local, *[k for k, v in LOCAL.items() if v == local]})] = local
# Atlantic cod does not occur here: a read tied between walleye pollock and G. morhua goes to
# pollock. Only at genus level, and only if pollock is a candidate - a read matching G. morhua
# alone is left at genus, since Pacific cod (G. macrocephalus, local) is Atlantic cod's closest
# relative.
GENUS_RULES[frozenset({"Gadus chalcogrammus", "Gadus finnmarchicus", "Gadus morhua"})] = "Gadus chalcogrammus"
# Sand lances (Ammodytidae) are assigned at family level only. The local species, Pacific sand
# lance (Ammodytes personatus), has few or no mitogenomes in the reference database, so its reads
# match other ammodytids best (e.g. Hyperoplus lanceolatus, a NE Atlantic species). Added Oct 8.
FAMILY_ONLY = {"Ammodytes": "Ammodytidae", "Hyperoplus": "Ammodytidae"}

# Chosen from the aquarium calibration (Sep 29): Q20, >=2 kb aligned, 98% species, delta 0.5.
DEFAULTS = dict(min_q=20.0, min_aln=2000, delta=0.5, sp_id=98.0, gen_id=95.0,
                fam_id=90.0, bleed=0.005, fragment_split=9200, assay_label="")


def min_aln_from(cfg):
    """Minimum aligned bp: ASSIGN_MIN_ALN in the config (MiFish: 120), else the LR-PCR default."""
    v = (cfg.get("ASSIGN_MIN_ALN") or "").strip()
    return int(v) if v else DEFAULTS["min_aln"]


# ----------------------------------------------------------------------------- config
def read_config(path):
    """Evaluate config.sh with bash and return its variables as a dict."""
    keys = ["PROJECT_DIR", "REF_DIR", "READ_DIR", "BLAST_DIR", "RESULT_DIR", "LOG_DIR",
            "AQ_SAMPLES", "FRAGMENT_SPLIT", "BLAST_MIN_Q", "ASSAY_LABEL", "ASSIGN_MIN_ALN",
            "REGIONAL_SPECIES"]
    script = f"source {shlex.quote(path)} >/dev/null 2>&1; " + "; ".join(
        f'printf "%s\\t%s\\n" {k} "${{{k}}}"' for k in keys)
    out = subprocess.run(["bash", "-c", script], capture_output=True, text=True, check=True).stdout
    # Split on newlines only and keep empty values: an unset key (e.g. ASSIGN_MIN_ALN in
    # config.sh) prints "KEY<TAB>" and must not lose its tab to strip().
    cfg = {}
    for line in out.split("\n"):
        if "\t" in line:
            k, v = line.split("\t", 1)
            cfg[k] = v.rstrip("\r")
    return cfg


# ----------------------------------------------------------------------------- loading
def load_lineages(path):
    """taxid_lineage.tsv from taxonkit: taxid<TAB>k;p;c;o;f;g;s"""
    lin = pd.read_csv(path, sep="\t", header=None, names=["taxid", "lineage"], dtype=str)
    parts = lin["lineage"].fillna("").str.split(";", expand=True).reindex(columns=range(7))
    parts.columns = RANKS
    parts = parts.replace("", np.nan)
    parts.insert(0, "taxid", pd.to_numeric(lin["taxid"], errors="coerce"))
    return parts.dropna(subset=["taxid"]).astype({"taxid": "int64"}).drop_duplicates("taxid")


BLAST_COLS = ["qseqid", "sseqid", "staxids", "pident", "length", "mismatch", "gapopen", "qstart",
              "qend", "sstart", "send", "evalue", "bitscore", "qlen", "slen", "qcovs"]


def aggregate_blast(path, chunksize=2_000_000):
    """Sum HSPs per (read, reference): aligned bp and identity weighted by length."""
    parts = []
    for ch in pd.read_csv(path, sep="\t", header=None, names=BLAST_COLS,
                          usecols=["qseqid", "sseqid", "staxids", "pident", "length"],
                          dtype={"qseqid": str, "sseqid": str, "staxids": str},
                          chunksize=chunksize):
        ch["taxid"] = pd.to_numeric(ch["staxids"].str.split(";").str[0], errors="coerce")
        ch["matches"] = ch["pident"] * ch["length"] / 100.0
        parts.append(ch.groupby(["qseqid", "sseqid", "taxid"], as_index=False)[["length", "matches"]].sum())
    if not parts:
        return pd.DataFrame(columns=["read", "ref", "taxid", "aln", "pid"])
    agg = pd.concat(parts).groupby(["qseqid", "sseqid", "taxid"], as_index=False)[["length", "matches"]].sum()
    agg["pid"] = 100.0 * agg["matches"] / agg["length"]
    return agg.rename(columns={"qseqid": "read", "sseqid": "ref", "length": "aln"})[
        ["read", "ref", "taxid", "aln", "pid"]]


def load_sample(sample, blast_dir, read_dir, cache_dir=None):
    """Aggregated hits + per-read info for one sample (cached after first load)."""
    ipath = os.path.join(read_dir, f"{sample}.readinfo.tsv")
    empty_hits = pd.DataFrame({"read": pd.Series(dtype=str), "ref": pd.Series(dtype=str),
                               "taxid": pd.Series(dtype="Int64"), "aln": pd.Series(dtype=float),
                               "pid": pd.Series(dtype=float)})
    # A sample with no reads (e.g. an extraction blank or PCR negative) is kept as zero reads:
    # BLAST is skipped for it, so there are no results to load.
    if not os.path.exists(ipath) or os.path.getsize(ipath) == 0:
        info = pd.DataFrame({"read": pd.Series(dtype=str), "qlen": pd.Series(dtype=int),
                             "meanq": pd.Series(dtype=float)})
        empty_hits.insert(0, "sample", sample); info.insert(0, "sample", sample)
        print(f"  {sample}: 0 reads (no BLAST needed)", file=sys.stderr)
        return empty_hits, info
    bpath = os.path.join(blast_dir, f"{sample}.blast.tsv.gz")
    if not os.path.exists(bpath):
        bpath = bpath[:-3]
    if not os.path.exists(bpath):
        raise FileNotFoundError(
            f"No BLAST results for {sample} in {blast_dir}. Run 03_blast_reads.sh for it first "
            f"(or, for the quick pilot, add --subsample N to both 03 and 05).")
    cache = os.path.join(cache_dir, f"{sample}.hits.tsv.gz") if cache_dir else None
    # Reuse the cached summary only if it is newer than the BLAST results it came from.
    if cache and os.path.exists(cache) and os.path.getmtime(cache) >= os.path.getmtime(bpath):
        hits = pd.read_csv(cache, sep="\t", dtype={"read": str, "ref": str})
    else:
        hits = aggregate_blast(bpath)
        if cache:
            os.makedirs(cache_dir, exist_ok=True)
            hits.to_csv(cache, sep="\t", index=False)
    info = pd.read_csv(ipath, sep="\t", header=None,
                       names=["read", "qlen", "meanq"], dtype={"read": str}, usecols=[0, 1, 2])
    hits.insert(0, "sample", sample)
    info.insert(0, "sample", sample)
    return hits, info


def list_samples(read_dir):
    return sorted(os.path.basename(p)[:-len(".readinfo.tsv")]
                  for p in glob.glob(os.path.join(read_dir, "*.readinfo.tsv")))


def load_all(samples, blast_dir, read_dir, cache_dir=None):
    H, I = zip(*(load_sample(s, blast_dir, read_dir, cache_dir) for s in samples))
    hits = pd.concat(H, ignore_index=True)
    hits["taxid"] = pd.to_numeric(hits["taxid"], errors="coerce").astype("Int64")
    return hits, pd.concat(I, ignore_index=True)


# ----------------------------------------------------------------------------- assignment
def load_regional(path):
    """One species name per line; blank lines and lines starting with # are ignored."""
    with open(path) as f:
        return {l.split("#")[0].strip() for l in f if l.split("#")[0].strip()}


def regional_path(cfg):
    """REGIONAL_SPECIES in the config, else lists/regional_species.txt in the project or script folder."""
    v = (cfg.get("REGIONAL_SPECIES") or "").strip()
    if v:
        return v
    for d in (cfg.get("PROJECT_DIR") or "", os.path.dirname(os.path.abspath(__file__))):
        f = os.path.join(d, "lists", "regional_species.txt")
        if d and os.path.exists(f):
            return f
    return None


def assign(hits, info, lineages, p, regional=None):
    """Return one row per read (passing min_q) with its assignment.

    regional: optional set of species names. A read assigned to a species that is not in the
    set is reported at genus level instead (e.g. NW Pacific Sebastes, Atlantic clupeids, human).
    """
    p = {**DEFAULTS, **p}
    # Reads without a quality value (FASTA input, e.g. the MiFish reads) are kept: their
    # quality filtering was done upstream, so the Q cut-off cannot be applied here.
    reads = info[(info["meanq"] >= p["min_q"]) | info["meanq"].isna()].copy()
    if p.get("assay_label"):            # single-amplicon assay (MiFish): one label for all reads
        reads["fragment"] = p["assay_label"]
    else:
        reads["fragment"] = np.where(reads["qlen"] < p["fragment_split"], "8.0kb", "9.5kb")
    key = ["sample", "read"]

    h = hits[hits["aln"] >= p["min_aln"]].merge(reads[key], on=key)
    if h.empty:
        reads[["rank", "taxon", "best_pid", "best_aln", "n_cand"]] = [None, None, np.nan, np.nan, 0]
        reads["status"] = "no_hit"
        return reads
    h["best_pid"] = h.groupby(key)["pid"].transform("max")
    cand = h[h["pid"] >= h["best_pid"] - p["delta"]].merge(lineages, on="taxid", how="left")

    g = cand.groupby(key, sort=False)
    size = g.size()
    # A rank agrees if every candidate has the same value there. A rank that is blank for
    # ALL candidates (e.g. "kingdom", which recent NCBI taxonomy leaves empty for fish)
    # is skipped rather than treated as a disagreement.
    agree = {}
    for r in RANKS:
        cnt, nun = g[r].count(), g[r].nunique()
        agree[r] = ((nun == 1) & (cnt == size)) | (cnt == 0)
    agree = pd.DataFrame(agree)
    depth = agree.cumprod(axis=1).sum(axis=1)                 # number of ranks agreed from the top
    summ = pd.DataFrame({"best_pid": g["best_pid"].first(),
                         "best_aln": g["aln"].max(),
                         "n_cand": g.size(),
                         "species_set": g["species"].agg(lambda s: frozenset(s.dropna())),
                         "depth": depth})
    first = g[RANKS].first()

    cap = np.select([summ["best_pid"] >= p["sp_id"], summ["best_pid"] >= p["gen_id"],
                     summ["best_pid"] >= p["fam_id"]], [7, 6, 5], default=0)
    summ["depth"] = np.minimum(summ["depth"], cap)
    # Step back to the deepest agreed rank that actually has a name.
    fv = first.to_numpy(dtype=object)
    depths = []
    for i, d in enumerate(summ["depth"].to_numpy()):
        d = int(d)
        while d > 0 and (fv[i, d - 1] is None or (isinstance(fv[i, d - 1], float) and np.isnan(fv[i, d - 1]))):
            d -= 1
        depths.append(d)
    summ["depth"] = depths
    summ["rank"] = summ["depth"].map(lambda d: RANKS[d - 1] if d > 0 else None)
    summ["taxon"] = [fv[i, d - 1] if d > 0 else None for i, d in enumerate(summ["depth"])]

    # Manuscript rules for flatfish congeners / synonyms (only where species is allowed).
    sp_ok = summ["best_pid"] >= p["sp_id"]
    for i in np.flatnonzero((summ["rank"] == "genus").to_numpy() & sp_ok.to_numpy()):
        s = summ["species_set"].iat[i]
        for members, name in GENUS_RULES.items():
            if s and s <= members and name in s:
                summ.iat[i, summ.columns.get_loc("rank")] = "species"
                summ.iat[i, summ.columns.get_loc("taxon")] = name
    is_sp = summ["rank"] == "species"
    summ.loc[is_sp, "taxon"] = summ.loc[is_sp, "taxon"].replace(SPECIES_REMAP)
    # Genera assigned at family level only (see FAMILY_ONLY).
    gen = summ["taxon"].astype(str).str.split().str[0]
    fo = summ["rank"].isin(["species", "genus"]) & gen.isin(FAMILY_ONLY)
    summ.loc[fo, "taxon"] = gen[fo].map(FAMILY_ONLY)
    summ.loc[fo, "rank"] = "family"
    # Regional species list (Oct 8): species-level calls to species not recorded in the region
    # are reported at genus level. Applied after the local-species rules above.
    if regional:
        nr = (summ["rank"] == "species") & ~summ["taxon"].isin(regional)
        summ["not_regional"] = np.where(nr, summ["taxon"], "")
        summ.loc[nr, "taxon"] = summ.loc[nr, "taxon"].astype(str).str.split().str[0]
        summ.loc[nr, "rank"] = "genus"

    # Keep the competing species for reads that stop above species level (diagnostics).
    summ["cand_species"] = [";".join(sorted(x)) if r in ("genus", "family") else ""
                            for x, r in zip(summ["species_set"], summ["rank"])]
    out = reads.merge(summ.drop(columns=["species_set", "depth"]).reset_index(), on=key, how="left")
    out["status"] = np.where(out["n_cand"].isna(), "no_hit",
                             np.where(out["rank"].isna(), "unassigned", "assigned"))
    out["n_cand"] = out["n_cand"].fillna(0).astype(int)
    return out


def count_taxa(assigned, bleed=DEFAULTS["bleed"]):
    """Reads per sample x fragment x taxon, with the barcode-bleed filter flagged."""
    a = assigned[assigned["status"] == "assigned"]
    c = a.groupby(["sample", "fragment", "rank", "taxon"], as_index=False).size().rename(
        columns={"size": "reads"})
    c["taxon_total"] = c.groupby(["fragment", "taxon"])["reads"].transform("sum")
    c["share_of_taxon"] = c["reads"] / c["taxon_total"]
    c["kept"] = c["share_of_taxon"] >= bleed
    return c.sort_values(["fragment", "sample", "reads"], ascending=[True, True, False])


def summarise(assigned, counts):
    s = assigned.groupby(["sample", "fragment"]).agg(
        reads=("read", "size"),
        no_hit=("status", lambda x: (x == "no_hit").sum()),
        unassigned=("status", lambda x: (x == "unassigned").sum()),
        to_species=("rank", lambda x: (x == "species").sum()),
        to_genus=("rank", lambda x: (x == "genus").sum()),
        to_family_or_above=("rank", lambda x: x.isin(RANKS[:5]).sum()),
        median_best_pid=("best_pid", "median"),
    ).reset_index()
    kept = counts[counts["kept"] & (counts["rank"] == "species")]
    n = kept.groupby(["sample", "fragment"]).size().rename("species_detected").reset_index()
    return s.merge(n, how="left").fillna({"species_detected": 0}).astype({"species_detected": int})


# ----------------------------------------------------------------------------- CLI
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", default=os.path.join(os.path.dirname(os.path.abspath(__file__)), "config.sh"))
    ap.add_argument("--samples", nargs="*", help="default: every sample in the reads folder")
    ap.add_argument("--out", help="output folder (default: RESULT_DIR/assignment)")
    ap.add_argument("--min-q", type=float, default=DEFAULTS["min_q"])
    ap.add_argument("--min-aln", type=int, default=None,
                    help="default: ASSIGN_MIN_ALN in the config, else 2000 (LR-PCR)")
    ap.add_argument("--delta", type=float, default=DEFAULTS["delta"])
    ap.add_argument("--sp-id", type=float, default=DEFAULTS["sp_id"])
    ap.add_argument("--gen-id", type=float, default=DEFAULTS["gen_id"])
    ap.add_argument("--fam-id", type=float, default=DEFAULTS["fam_id"])
    ap.add_argument("--bleed", type=float, default=DEFAULTS["bleed"])
    ap.add_argument("--regional", help="regional species list (default: REGIONAL_SPECIES in the config, "
                    "else lists/regional_species.txt)")
    ap.add_argument("--no-regional", action="store_true", help="do not apply the regional species list")
    a = ap.parse_args()

    cfg = read_config(a.config)
    out = a.out or os.path.join(cfg["RESULT_DIR"], "assignment")
    os.makedirs(out, exist_ok=True)
    samples = a.samples or list_samples(cfg["READ_DIR"])
    assay = (cfg.get("ASSAY_LABEL") or "").strip()
    min_aln = a.min_aln if a.min_aln is not None else min_aln_from(cfg)
    if assay:
        # MiFish: every sample file in the reads folder is used (names are not barcodeNN).
        # Keep one file per sample there; a pooled file would distort the barcode-bleed filter.
        meta_path = None
    # Only real barcodes: a pooled file (e.g. all samples concatenated) left in the reads folder
    # would double every taxon's total and break the barcode-bleed filter.
    else:
        meta_path = os.path.join(cfg["PROJECT_DIR"], "lists", "sample_metadata.csv")
    if assay:
        keep = list(samples)
    elif os.path.exists(meta_path):
        known = set(pd.read_csv(meta_path)["demulti"].astype(str))
        keep = [x for x in samples if (re.search(r"barcode\d+", x) or [None])[0] in known]
    else:
        keep = [x for x in samples if re.search(r"barcode\d+", x)]
    for x in sorted(set(samples) - set(keep)):
        print(f"  skipping {x}: not a barcode in lists/sample_metadata.csv", file=sys.stderr)
    samples = keep
    if not samples:
        sys.exit("No samples to process: LR-PCR sample names must contain a barcode listed in "
                 "lists/sample_metadata.csv (MiFish: check READ_DIR in config_mifish.sh).")
    params = dict(min_q=a.min_q, min_aln=min_aln, delta=a.delta, sp_id=a.sp_id, gen_id=a.gen_id,
                  fam_id=a.fam_id, bleed=a.bleed, fragment_split=int(cfg.get("FRAGMENT_SPLIT") or 9200),
                  assay_label=assay)

    lineages = load_lineages(os.path.join(cfg["REF_DIR"], "taxid_lineage.tsv"))
    hits, info = load_all(samples, cfg["BLAST_DIR"], cfg["READ_DIR"], os.path.join(cfg["RESULT_DIR"], "cache"))
    reg_file = None if a.no_regional else (a.regional or regional_path(cfg))
    regional = load_regional(reg_file) if reg_file else None
    params["regional_species"] = reg_file or ""
    print(f"Regional species list: {reg_file or 'none (all species-level calls kept)'}"
          + (f" ({len(regional)} species)" if regional else ""), file=sys.stderr)
    assigned = assign(hits, info, lineages, params, regional)
    if regional and "not_regional" in assigned:
        # Report what the regional rule moved to genus, so the list can be checked.
        moved = assigned[assigned["not_regional"].fillna("").ne("")]
        if len(moved):
            print("Species not on the regional list, reported at genus level (reads):", file=sys.stderr)
            for sp_name, n in moved["not_regional"].value_counts().items():
                print(f"  {sp_name}\t{n}", file=sys.stderr)
    counts = count_taxa(assigned, params["bleed"])
    summary = summarise(assigned, counts)
    # Keep samples with no reads (controls) in the summary as explicit zeros.
    missing = [x for x in samples if x not in set(summary["sample"])]
    if missing:
        z = pd.DataFrame({"sample": missing, "fragment": "none", "reads": 0, "no_hit": 0, "unassigned": 0,
                          "to_species": 0, "to_genus": 0, "to_family_or_above": 0,
                          "median_best_pid": np.nan, "species_detected": 0})
        summary = pd.concat([summary, z], ignore_index=True).sort_values(["sample", "fragment"])

    assigned.to_csv(os.path.join(out, "read_assignments.tsv.gz"), sep="\t", index=False)
    counts.to_csv(os.path.join(out, "taxon_counts.tsv"), sep="\t", index=False)
    summary.to_csv(os.path.join(out, "sample_summary.tsv"), sep="\t", index=False)
    sp = counts[counts["kept"] & (counts["rank"] == "species")]
    sp.pivot_table(index=["fragment", "taxon"], columns="sample", values="reads", fill_value=0) \
      .to_csv(os.path.join(out, "species_by_sample.tsv"), sep="\t")
    # Which species did genus-level reads fall between? (e.g. Merluccius, Gadus, Sebastes)
    gl = assigned[assigned["rank"] == "genus"]
    if len(gl):
        (gl.groupby(["taxon", "cand_species"]).size().rename("reads").reset_index()
           .sort_values(["taxon", "reads"], ascending=[True, False])
           .to_csv(os.path.join(out, "genus_level_candidates.tsv"), sep="\t", index=False))
    json.dump(params, open(os.path.join(out, "params.json"), "w"), indent=2)
    print(summary.to_string(index=False))
    print(f"\nWrote results to {out}")


if __name__ == "__main__":
    main()
