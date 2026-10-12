#!/usr/bin/env python3
"""
Download sequences + metadata for a list of NCBI nucleotide GIs, straight from the
E-utilities web service (no EDirect needed). Called by 01_build_refdb.sh.

Writes, per chunk of 500 GIs:  <outdir>/c_NNNNN.fa   (FASTA)
                               <outdir>/c_NNNNN.tsv  (accession.version, taxid, length, title)
Finished chunks are skipped, so an interrupted run resumes where it stopped.
Respects NCBI's rate limit: 3 requests/s, or 10/s with NCBI_API_KEY set.

Usage: python3 fetch_ncbi.py all_gis.txt chunks/
"""
import json, os, ssl, sys, time, threading, urllib.error, urllib.parse, urllib.request
from concurrent.futures import ThreadPoolExecutor, as_completed

BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
CHUNK = 500
KEY = os.environ.get("NCBI_API_KEY", "").strip()
RATE = 9 if KEY else 2.5            # requests per second, a little under NCBI's limit
WORKERS = 4 if KEY else 2

# NCBI closes connections without a TLS close_notify; OpenSSL 3 (as shipped by conda)
# treats that as an error ("unexpected eof while reading"). Accept it, and instead
# verify each chunk is complete (FASTA record count == esummary record count).
_CTX = ssl.create_default_context()
_CTX.options |= getattr(ssl, "OP_IGNORE_UNEXPECTED_EOF", 0)

_lock = threading.Lock()
_last = [0.0]


def _throttle():
    with _lock:
        wait = _last[0] + 1.0 / RATE - time.time()
        if wait > 0:
            time.sleep(wait)
        _last[0] = time.time()


def post(endpoint, params, tries=6):
    params = dict(params, tool="lrpcr_blast")
    if KEY:
        params["api_key"] = KEY
    data = urllib.parse.urlencode(params).encode()
    for t in range(tries):
        _throttle()
        try:
            with urllib.request.urlopen(urllib.request.Request(BASE + endpoint, data=data), timeout=300, context=_CTX) as r:
                return r.read().decode("utf-8", errors="replace")
        except (urllib.error.URLError, OSError) as e:  # includes SSL and timeout errors
            code = getattr(e, "code", None)
            pause = 30 if code == 429 else 10 * (t + 1)
            print(f"    {endpoint} error ({code or e}); retry in {pause}s", file=sys.stderr, flush=True)
            time.sleep(pause)
    raise RuntimeError(f"{endpoint} failed after {tries} attempts")


def fasta_lengths(text):
    """accession -> sequence length for every record in a FASTA string."""
    out, acc, n = {}, None, 0
    for line in text.splitlines():
        if line.startswith(">"):
            if acc is not None:
                out[acc] = n
            acc, n = line[1:].split(None, 1)[0] if len(line) > 1 else "", 0
        else:
            n += len(line.strip())
    if acc is not None:
        out[acc] = n
    return out


def keep_records(text, keep):
    """Return only the FASTA records whose accession is in `keep`."""
    out, on = [], False
    for line in text.splitlines(keepends=True):
        if line.startswith(">"):
            on = (line[1:].split(None, 1)[0] if len(line) > 1 else "") in keep
        if on:
            out.append(line)
    return "".join(out)


def fetch_chunk(i, gis, outdir):
    fa, tsv = os.path.join(outdir, f"c_{i:05d}.fa"), os.path.join(outdir, f"c_{i:05d}.tsv")
    if os.path.exists(fa) and os.path.exists(tsv):     # written only once a chunk is complete
        return i, "skip", 0
    ids = ",".join(gis)
    prev_bad = None
    for attempt in range(5):
        try:
            summ = json.loads(post("esummary.fcgi", {"db": "nuccore", "id": ids, "retmode": "json"}))
        except json.JSONDecodeError:
            continue                                   # truncated reply: fetch again
        res = summ.get("result", {})
        rows, expected = {}, {}
        for uid in res.get("uids", []):
            d = res.get(uid, {})
            if "error" in d or not d.get("accessionversion"):
                continue
            title = str(d.get("title", "")).replace("\t", " ")
            rows[d["accessionversion"]] = f"{d['accessionversion']}\t{d.get('taxid', '')}\t{d.get('slen', '')}\t{title}"
            expected[d["accessionversion"]] = int(d.get("slen") or 0)
        fasta = post("efetch.fcgi", {"db": "nuccore", "id": ids, "rettype": "fasta", "retmode": "text"})
        got = fasta_lengths(fasta)
        # Complete = every record NCBI lists is present at its full length.
        bad = {a: got.get(a) for a, n in expected.items() if got.get(a) != n}
        if expected and not bad:
            break
        # A cut-off download differs from try to try. The same records missing twice in a
        # row means NCBI itself does not serve them (removed/replaced/suppressed records):
        # accept the chunk, drop those records, and log them.
        if bad and bad == prev_bad:
            with _lock, open(os.path.join(outdir, "unavailable_records.tsv"), "a") as f:
                for a, n in bad.items():
                    f.write(f"{i}\t{a}\t{expected[a]}\t{'missing' if n is None else n}\n")
            print(f"    chunk {i}: {len(bad)} record(s) not served by NCBI; logged and skipped",
                  file=sys.stderr, flush=True)
            break
        prev_bad = bad
        print(f"    chunk {i}: {len(bad)} of {len(expected)} records missing or truncated; "
              f"refetching", file=sys.stderr, flush=True)
        time.sleep(5)
    else:
        raise RuntimeError(f"chunk {i}: still incomplete after 5 attempts")
    # Keep only complete records in both files so sequences and metadata match.
    keep = {a for a in expected if a not in bad}
    fasta = keep_records(fasta, keep)
    with open(fa + ".tmp", "w") as f:
        f.write(fasta)
    with open(tsv + ".tmp", "w") as f:
        f.write("\n".join(rows[a] for a in rows if a in keep) + "\n")
    os.replace(fa + ".tmp", fa)
    os.replace(tsv + ".tmp", tsv)
    return i, "ok", fasta.count(">")


def main():
    gi_file, outdir = sys.argv[1], sys.argv[2]
    os.makedirs(outdir, exist_ok=True)
    gis = [l.strip() for l in open(gi_file) if l.strip().isdigit()]
    chunks = [gis[k:k + CHUNK] for k in range(0, len(gis), CHUNK)]
    print(f"{len(gis)} GIs in {len(chunks)} chunks; {WORKERS} workers at <= {RATE}/s "
          f"({'with' if KEY else 'without'} API key)", flush=True)
    done = failed = 0
    t0 = time.time()
    with ThreadPoolExecutor(WORKERS) as ex:
        futs = {ex.submit(fetch_chunk, i, c, outdir): i for i, c in enumerate(chunks)}
        for fut in as_completed(futs):
            try:
                i, status, n = fut.result()
                done += 1
            except Exception as e:
                failed += 1
                print(f"  chunk {futs[fut]} FAILED: {e}", file=sys.stderr, flush=True)
                continue
            if done % 20 == 0 or done == len(chunks):
                el = time.time() - t0
                print(f"  {done}/{len(chunks)} chunks ({el / 60:.0f} min elapsed)", flush=True)
    # Record the GIs of any chunk that is still not on disk (failed now or earlier).
    missing = [g for i, c in enumerate(chunks)
               if not os.path.exists(os.path.join(outdir, f"c_{i:05d}.tsv")) for g in c]
    miss_path = os.path.join(outdir, "missing_gis.txt")
    with open(miss_path, "w") as f:
        f.write("".join(g + "\n" for g in missing))
    if missing:
        print(f"{len(missing)} of {len(gis)} GIs ({100 * len(missing) / len(gis):.2f}%) not downloaded; "
              f"listed in {miss_path}", file=sys.stderr, flush=True)
        sys.exit(2)


if __name__ == "__main__":
    main()
