#!/usr/bin/env python3
"""
Step 6 - reference coverage for the aquarium species (R1 major 2; R2 8).

For each species in lists/aquarium_species.tsv, counts records in the reference
database: all records, complete mitogenomes, and other records at least --min-len long
(default 2,000 bp, the minimum aligned length used for LR-PCR reads). A tank species with no
complete mitogenome and no record that long cannot reach the minimum alignment, which
separates "missing reference" from "not amplified".

Input: a metadata table of the BLAST database, one row per record, tab-separated:
    accession, taxid, length, scientific name, title
Make it from any BLAST database with blastdbcmd (as run_06_reference_coverage.slurm does):
    blastdbcmd -db vert_mtDNA -entry all -outfmt "%a	%T	%l	%S	%t" > ref_metadata.tsv
The older 4-column ref_metadata.tsv from 01_build_refdb.sh (no scientific name) also works
if --config points to a REF_DIR that has it and taxid_lineage.tsv.

The species of each record comes from the scientific name column (BLAST taxonomy database),
or, when that is missing ("N/A"), from the first two words of the record title.

Usage:  python3 06_reference_coverage.py --metadata ref_metadata.tsv --out reference_coverage.tsv
        python3 06_reference_coverage.py            # old layout: REF_DIR/ref_metadata.tsv via config.sh
Output: reference_coverage.tsv (default: RESULT_DIR/reference_coverage.tsv)
"""
import argparse, os, re, sys
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))


def species_from_title(t):
    t = re.sub(r"^(PREDICTED|UNVERIFIED|TPA|TSA|MAG):\s*", "", str(t), flags=re.I)
    w = t.split()
    return f"{w[0]} {w[1]}" if len(w) >= 2 and w[0][:1].isupper() and w[1][:1].islower() else None


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--metadata", help="blastdbcmd table: accession, taxid, length, sci name, title")
    ap.add_argument("--config", default=os.path.join(HERE, "config.sh"))
    ap.add_argument("--species", default=os.path.join(HERE, "lists", "aquarium_species.tsv"))
    ap.add_argument("--out", help="output file (default: RESULT_DIR/reference_coverage.tsv)")
    ap.add_argument("--min-len", type=int, default=2000, help="length cut-off for other long records (bp)")
    a = ap.parse_args()

    if a.metadata:
        meta = pd.read_csv(a.metadata, sep="\t", header=None, names=["acc", "taxid", "slen", "sciname", "title"],
                           dtype=str, quoting=3, on_bad_lines="skip")
        meta["slen"] = pd.to_numeric(meta["slen"], errors="coerce").fillna(0)
        # species = first two words of the scientific name (drops subspecies), else from the title.
        # Without the BLAST taxonomy files every scientific name is "N/A".
        def two_words(x):
            w = str(x).split() if isinstance(x, str) else []
            return " ".join(w[:2]) if len(w) >= 2 and x.strip() not in ("N/A", "NA") else None
        meta["species"] = meta["sciname"].map(two_words)
        miss = meta["species"].isna()
        meta.loc[miss, "species"] = meta.loc[miss, "title"].map(species_from_title)
        print(f"{len(meta):,} records; species from BLAST taxonomy for {(~miss).sum():,}, "
              f"from titles for {miss.sum():,}", file=sys.stderr)
        out_path = a.out or "reference_coverage.tsv"
    else:
        sys.path.insert(0, HERE)
        import lrpcr_assign as la
        cfg = la.read_config(a.config)
        meta = pd.read_csv(os.path.join(cfg["REF_DIR"], "ref_metadata.tsv"), sep="\t", header=None,
                           names=["acc", "taxid", "slen", "title"], dtype={"acc": str, "title": str},
                           on_bad_lines="skip")
        meta["taxid"] = pd.to_numeric(meta["taxid"], errors="coerce")
        lin = la.load_lineages(os.path.join(cfg["REF_DIR"], "taxid_lineage.tsv"))[["taxid", "species"]]
        meta = meta.merge(lin, on="taxid", how="left")
        out_path = a.out or os.path.join(cfg["RESULT_DIR"], "reference_coverage.tsv")

    t = meta["title"].fillna("").str.lower()
    meta["mitogenome"] = t.str.contains("complete genome|complete mitochondrial genome") & (meta["slen"] >= 15000)
    meta["long_other"] = ~meta["mitogenome"] & (meta["slen"] >= a.min_len)

    sp = pd.read_csv(a.species, sep="\t")
    other = f"other_ge_{a.min_len // 1000}kb" if a.min_len % 1000 == 0 else f"other_ge_{a.min_len}bp"
    agg = meta.groupby("species").agg(records=("acc", "size"), mitogenomes=("mitogenome", "sum"),
                                      **{other: ("long_other", "sum")}).reset_index()
    out = sp.merge(agg, on="species", how="left").fillna(0)
    for c in ["records", "mitogenomes", other]:
        out[c] = out[c].astype(int)
    d = os.path.dirname(out_path)
    if d:
        os.makedirs(d, exist_ok=True)
    out.to_csv(out_path, sep="\t", index=False)
    print(out.to_string(index=False))
    n = out[out["tank_census"] > 0]
    print(f"\nTank species (census > 0): {len(n)}; with >=1 complete mitogenome: {(n['mitogenomes'] > 0).sum()}; "
          f"with no record >= {a.min_len:,} bp: {((n['mitogenomes'] + n[other]) == 0).sum()}")
    print(f"Wrote {out_path}")


if __name__ == "__main__":
    main()
