#!/usr/bin/env python3
r"""Fetch ISMIP7 AIS forcing from the Source Cooperative mirror over plain
HTTPS (anonymous, resumable), no AWS tooling needed.

    python antarctica/scripts/download_mirror.py PREFIX [PREFIX ...]
        [--root ISMIP7/AIS] [--include REGEX] [--dry-run] [--jobs 4]

PREFIX is product-relative, e.g. ``data/CESM2-WACCM/ctrl/ocean/tf/`` or
``data/OCX/ocean/main/``. Files land under ``--root`` with the mirror's
layout minus the leading ``data/`` (so ``data/CESM2-WACCM/ctrl/...`` becomes
``<root>/CESM2-WACCM/ctrl/...``), and versioned files are placed in a
``<version>/`` directory the way the Globus tree and this repository's
readers expect (the mirror itself keeps no version directories).

Partial files resume with a Range request; a file whose size matches the
mirror is skipped. Every transfer is checked against the size the listing
gave: a short one is deleted and reported FAILED rather than left on disk as
a truncated NetCDF. The mirror refuses Python's default User-Agent.
"""
import argparse
import concurrent.futures as cf
import os
import re
import sys
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET

ENDPOINT = "https://data.source.coop/ismip/ismip7-ais-forcing/"
PRODUCT = "ismip7-ais-forcing/"
HEADERS = {"User-Agent": "curl/8"}
VERSION = re.compile(r"[_-](v\d+(?:\.\d+)*)(?:_|\.nc$)")


def list_keys(prefix):
    keys, token = [], None
    while True:
        url = ENDPOINT + "?list-type=2&max-keys=1000&prefix=" + urllib.parse.quote(prefix)
        if token:
            url += "&continuation-token=" + urllib.parse.quote(token)
        with urllib.request.urlopen(urllib.request.Request(url, headers=HEADERS), timeout=60) as r:
            root = ET.fromstring(r.read())
        ns = {"s3": root.tag.split("}")[0].strip("{")}
        for k in root.findall("s3:Contents", ns):
            keys.append((k.find("s3:Key", ns).text[len(PRODUCT):], int(k.find("s3:Size", ns).text)))
        nxt = root.find("s3:NextContinuationToken", ns)
        if nxt is None:
            return keys
        token = nxt.text


def local_path(root, key):
    r"""``data/<esm>/<scenario>/<product>/<variable>/<file>`` ->
    ``<root>/<esm>/<scenario>/<product>/<variable>/<version>/<file>``."""
    rel = key[len("data/"):] if key.startswith("data/") else key
    d, f = os.path.split(rel)
    m = VERSION.search(f)
    if m and os.path.basename(d) != m.group(1):
        d = os.path.join(d, m.group(1))
    return os.path.join(root, d, f)


def fetch(key, size, dest):
    os.makedirs(os.path.dirname(dest), exist_ok=True)
    have = os.path.getsize(dest) if os.path.exists(dest) else 0
    if have == size:
        return "skip"
    if have > size:
        os.remove(dest); have = 0
    req = urllib.request.Request(ENDPOINT + urllib.parse.quote(key[len(PRODUCT):] if key.startswith(PRODUCT) else key), headers=dict(HEADERS))
    if have:
        req.add_header("Range", f"bytes={have}-")
    with urllib.request.urlopen(req, timeout=120) as r, open(dest, "ab" if have else "wb") as out:
        while True:
            chunk = r.read(1 << 22)
            if not chunk:
                break
            out.write(chunk)
    # An empty chunk from a dropped stream is indistinguishable from the end
    # of the body, so a truncated NetCDF would otherwise land on disk and be
    # reported as fetched; the readers only discover it years of forcing
    # later. Compare against the size the listing gave and fail loudly.
    got = os.path.getsize(dest)
    if got != size:
        os.remove(dest)
        raise IOError(f"short transfer: got {got} of {size} bytes; removed the partial file")
    return "resumed" if have else "fetched"


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("prefix", nargs="+")
    ap.add_argument("--root", default=os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))), "ISMIP7", "AIS"))
    ap.add_argument("--include", default=None, help="regex a key must match")
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--jobs", type=int, default=4)
    a = ap.parse_args()
    todo = []
    for p in a.prefix:
        for key, size in list_keys(p):
            if a.include and not re.search(a.include, key):
                continue
            todo.append((key, size, local_path(a.root, key)))
    total = sum(s for _, s, _ in todo)
    print(f"{len(todo)} files, {total / 1e9:.2f} GB -> {a.root}", flush=True)
    if a.dry_run:
        for key, size, dest in todo[:20]:
            print(f"  {size / 1e6:9.1f} MB  {key}  ->  {os.path.relpath(dest, a.root)}")
        if len(todo) > 20:
            print(f"  ... {len(todo) - 20} more")
        return 0
    done = 0
    with cf.ThreadPoolExecutor(max_workers=a.jobs) as ex:
        futs = {ex.submit(fetch, key, size, dest): (key, size) for key, size, dest in todo}
        for fut in cf.as_completed(futs):
            key, size = futs[fut]
            try:
                status = fut.result()
            except Exception as e:
                print(f"  FAILED {key}: {e}", flush=True); continue
            done += size
            print(f"  {status:8s} {size / 1e6:9.1f} MB  {key}   [{done / 1e9:.2f} of {total / 1e9:.2f} GB]", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
