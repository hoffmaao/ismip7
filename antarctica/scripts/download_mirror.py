#!/usr/bin/env python3
r"""Fetch ISMIP7 AIS data from the Source Cooperative mirror over plain
HTTPS (anonymous, resumable), no AWS tooling needed.

    python antarctica/scripts/download_mirror.py PREFIX [PREFIX ...]
        [--product ismip7-ais-forcing] [--root ISMIP7/AIS]
        [--include REGEX] [--dry-run] [--check] [--jobs 4]

``--product`` selects the mirror product. The default is the forcing; the
observations MIPkit lives in ``ismip7-ais-observations``.

PREFIX is product-relative, e.g. ``data/CESM2-WACCM/ctrl/ocean/tf/`` or
``data/OCX/ocean/main/``. Files land under ``--root`` where this repository's
readers look for them: the mirror's layout minus the leading ``data/``, under
the product's own root (the observations product lives in ``obs/``), with
versioned files placed in a ``<version>/`` directory the way the Globus tree
expects. The mirror itself keeps no version directories.

Partial files resume with a Range request. A complete file is skipped only
when it is still the object the mirror publishes: the focus groups replace
files in place, same name and same version (the shifted OCX ``dacabfdz`` of
discussion #45, the wrong MRI ``lake_properties`` of #41), and a size check
alone never sees that. So the ETag of every fetched object is kept in
``<root>/.mirror_manifest.json``, and a file whose ETag has moved is REPLACED:
fetched again beside the old copy and swapped in once complete.

A complete file the manifest has never seen (anything fetched before it
existed, or over Globus) carries no ETag to compare. If it landed after the
mirror's object was last written it is current and is adopted into the
manifest. If the mirror's object is newer it is OLDER, which proves nothing:
most of a Globus tree predates the mirror itself. OLDER files are listed and
left alone unless ``--older refetch`` (fetch them again) or ``--older adopt``
(vouch for them) says otherwise. ``--dry-run`` lists both kinds.

The totals line counts everything listed, which says nothing of how much is
left to move, so the plan is printed by verb under it, followed by every row
(``<ESM>/<scenario>/<product>/<variable>``, the audit's unit) that holds
anything other than ``skip`` or ``adopt``. ``--dry-run`` stops there: it
transfers nothing and writes nothing, the manifest included. ``--check`` is a
dry run whose exit status is 1 unless every file is ``skip`` or ``adopt``,
which is the file-level completeness gate ``audit_forcing_versions.py`` cannot
be (issue #41).

Every transfer is checked against the size the listing gave, so a truncated
NetCDF is never reported as fetched: a short one is left in place for the
next run to resume and reported FAILED, and one longer than the listing (a
server that ignored the Range header) is deleted and refetched. Either way
the exit status is non-zero. The mirror refuses Python's default User-Agent.
"""
import argparse
import concurrent.futures as cf
import datetime
import json
import os
import re
import sys
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET

MIRROR = "https://data.source.coop/ismip/"
DEFAULT_PRODUCT = "ismip7-ais-forcing"
HEADERS = {"User-Agent": "curl/8"}
VERSION = re.compile(r"[_-](v\d+(?:\.\d+)*)(?:_|\.nc$)")

# Where a product's tree sits under --root. The forcing is the root itself;
# the observations product is the obs/ subtree the readers search.
DEST_ROOT = {"ismip7-ais-observations/": "obs"}

# Destinations the Globus tree keeps flat, so no <version>/ directory is
# inserted. obs_dhdt._obs_kit_path lists obs/mipkit non-recursively and picks
# the newest MIPkit by the version in its filename.
FLAT_DESTS = ("obs/mipkit",)

MANIFEST_NAME = ".mirror_manifest.json"

# Every verdict ``plan`` returns, and the two that leave nothing to transfer.
VERBS = ("skip", "adopt", "fetch", "resume", "REPLACED", "OLDER")
CLEAN = ("skip", "adopt")


def list_keys(prefix, endpoint, product):
    r"""``[(key, size, etag, last_modified)]`` under ``prefix``, paginated.
    ``last_modified`` is seconds since the epoch."""
    keys, token = [], None
    while True:
        url = endpoint + "?list-type=2&max-keys=1000&prefix=" + urllib.parse.quote(prefix)
        if token:
            url += "&continuation-token=" + urllib.parse.quote(token)
        with urllib.request.urlopen(urllib.request.Request(url, headers=HEADERS), timeout=60) as r:
            root = ET.fromstring(r.read())
        ns = {"s3": root.tag.split("}")[0].strip("{")}
        for k in root.findall("s3:Contents", ns):
            key = k.find("s3:Key", ns).text
            if key.startswith(product):
                key = key[len(product):]
            etag = k.find("s3:ETag", ns)
            stamp = k.find("s3:LastModified", ns)
            keys.append((key, int(k.find("s3:Size", ns).text),
                         etag.text.strip('"') if etag is not None else "",
                         _epoch(stamp.text) if stamp is not None else 0.0))
        nxt = root.find("s3:NextContinuationToken", ns)
        if nxt is None:
            return keys
        token = nxt.text


def _epoch(stamp):
    r"""An S3 ``LastModified`` (``2026-08-29T12:08:52.000Z``) as epoch seconds."""
    return datetime.datetime.fromisoformat(stamp.replace("Z", "+00:00")).timestamp()


def load_manifest(root):
    r"""``{product + key: {size, etag, last_modified}}``, empty when absent."""
    try:
        with open(os.path.join(root, MANIFEST_NAME)) as f:
            return json.load(f)
    except (OSError, ValueError):
        return {}


def save_manifest(root, manifest):
    os.makedirs(root, exist_ok=True)
    tmp = os.path.join(root, MANIFEST_NAME + ".tmp")
    with open(tmp, "w") as f:
        json.dump(manifest, f, indent=1, sort_keys=True)
    os.replace(tmp, os.path.join(root, MANIFEST_NAME))


def plan(size, etag, last_modified, dest, entry):
    r"""What a file on disk needs: ``fetch``, ``resume``, ``skip``, ``adopt``,
    ``REPLACED`` or ``OLDER``. ``entry`` is the manifest's record, or None.

    ``adopt`` is a complete file the manifest has never seen and the mirror has
    not touched since it landed: it is current, and only needs recording.
    ``OLDER`` is such a file that landed before the mirror's object was
    written, which may or may not be the same bytes.
    """
    have = os.path.getsize(dest) if os.path.exists(dest) else 0
    if have == 0:
        return "fetch"
    if have != size:
        # Once whole and recorded, now the wrong length: the object changed.
        if entry:
            return "REPLACED"
        # A prefix of an object that has since been replaced is a prefix of
        # nothing, so only one the mirror has not touched since is resumed.
        if have < size and os.path.getmtime(dest) >= last_modified:
            return "resume"
        return "fetch"
    if entry:
        return "skip" if not etag or entry.get("etag") in ("", etag) else "REPLACED"
    return "OLDER" if os.path.getmtime(dest) < last_modified else "adopt"


def row_of(key):
    r"""The audit's row a key sits under: ``<ESM>/<scenario>/<product>/<variable>``,
    or as much of it as the key has (a flat fracture directory stops at the
    product)."""
    parts = key.split("/")            # data, esm, scenario, [product], [variable], ..., file
    return "/".join(parts[1:min(len(parts) - 1, 5)])


def summarize(todo):
    r"""``({verb: [files, bytes]}, {(verb, row): [files, bytes]})`` over a plan
    of ``(key, size, etag, stamp, dest, action)``. Every verb is in the first,
    so a zero prints as a zero. The second holds only what is left to look at:
    the rows with anything other than ``skip`` or ``adopt`` under them."""
    verbs = {v: [0, 0] for v in VERBS}
    rows = {}
    for key, size, _, _, _, action in todo:
        hits = [verbs[action]]
        if action not in CLEAN:
            hits.append(rows.setdefault((action, row_of(key)), [0, 0]))
        for h in hits:
            h[0] += 1
            h[1] += size
    return verbs, rows


def local_path(root, key, dest_root=""):
    r"""``data/<esm>/<scenario>/<product>/<variable>/<file>`` ->
    ``<root>/<esm>/<scenario>/<product>/<variable>/<version>/<file>``, and for
    the observations product ``data/mipkit/<file>`` -> ``<root>/obs/mipkit/<file>``."""
    rel = key[len("data/"):] if key.startswith("data/") else key
    if dest_root:
        rel = os.path.join(dest_root, rel)
    d, f = os.path.split(rel)
    m = VERSION.search(f)
    if m and os.path.basename(d) != m.group(1) and d not in FLAT_DESTS:
        d = os.path.join(d, m.group(1))
    return os.path.join(root, d, f)


def fetch(key, size, dest, endpoint, action="fetch"):
    r"""Carry out ``action`` (from ``plan``) and return the status to report."""
    if action in ("skip", "adopt", "OLDER"):
        return action
    os.makedirs(os.path.dirname(dest), exist_ok=True)
    final = dest
    if action == "REPLACED":
        # The copy on disk is complete but stale. It stays where the readers
        # find it until its replacement is whole, so a dropped transfer leaves
        # a run with yesterday's file rather than a truncated one.
        dest = final + ".new"
    if action != "resume" and os.path.exists(dest):
        os.remove(dest)
    have = os.path.getsize(dest) if os.path.exists(dest) else 0
    req = urllib.request.Request(endpoint + urllib.parse.quote(key), headers=dict(HEADERS))
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
    # later. Compare against the size the listing gave.
    got = os.path.getsize(dest)
    if got < size and dest != final:
        # A replacement is fetched whole or not at all: the object may move
        # again before the next run, so its prefix is not worth keeping.
        os.remove(dest)
        raise IOError(f"short transfer: {got} of {size} bytes; the stale copy stays in place")
    if got < size:
        # A dropped stream: the bytes on disk are a valid prefix, so they stay
        # for the next run's Range request. The transfer is still a failure,
        # which is what carries it into the exit status.
        raise IOError(f"short transfer: {got} of {size} bytes, kept for the next resume")
    if got > size:
        # The server ignored the Range header and appended a second copy of
        # the body, so what is on disk is unusable and cannot be resumed.
        os.remove(dest)
        raise IOError(f"over-long transfer: {got} of {size} bytes; removed, it will be refetched")
    if dest != final:
        os.replace(dest, final)
        return "REPLACED"
    return "resumed" if have else "fetched"


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("prefix", nargs="+")
    ap.add_argument("--product", default=DEFAULT_PRODUCT,
                    help="mirror product (default %(default)s; the "
                         "observations MIPkit is ismip7-ais-observations)")
    ap.add_argument("--root", default=os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))), "ISMIP7", "AIS"))
    ap.add_argument("--include", default=None, help="regex a key must match")
    ap.add_argument("--older", choices=("report", "refetch", "adopt"), default="report",
                    help="what to do with a complete file that predates the manifest "
                         "and is older than the mirror's object (default %(default)s)")
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--check", action="store_true",
                    help="a dry run whose exit status is 1 unless every file is skip or adopt")
    ap.add_argument("--jobs", type=int, default=4)
    a = ap.parse_args()
    product = a.product.strip("/") + "/"
    endpoint = MIRROR + product
    dest_root = DEST_ROOT.get(product, "")
    manifest = load_manifest(a.root)
    todo = []
    for p in a.prefix:
        for key, size, etag, stamp in list_keys(p, endpoint, product):
            if a.include and not re.search(a.include, key):
                continue
            dest = local_path(a.root, key, dest_root)
            action = plan(size, etag, stamp, dest, manifest.get(product + key))
            if action == "OLDER" and a.older != "report":
                action = {"refetch": "REPLACED", "adopt": "adopt"}[a.older]
            todo.append((key, size, etag, stamp, dest, action))
    total = sum(t[1] for t in todo)
    print(f"{len(todo)} files, {total / 1e9:.2f} GB  {product} -> {a.root}", flush=True)
    verbs, rows = summarize(todo)
    print("by verb: " + ", ".join(f"{n} {v}" + (f" ({b / 1e9:.2f} GB)" if n else "")
                                  for v, (n, b) in verbs.items()), flush=True)
    for (verb, row), (n, b) in sorted(rows.items()):
        print(f"  {verb:8s} {n:6d} files {b / 1e9:9.2f} GB  {row}", flush=True)
    for verdict, what in (("REPLACED", "replaced on the mirror since they were fetched "
                                       "(same name, new content)"),
                          ("OLDER", "older than the mirror's object and unknown to the "
                                    "manifest; left alone, see --older")):
        hits = [t for t in todo if t[5] == verdict]
        if hits:
            print(f"{len(hits)} {what}:", flush=True)
            for key, *_ in hits:
                print(f"  {verdict:8s} {key}", flush=True)
    if a.dry_run or a.check:
        for key, size, _, _, dest, action in todo[:20]:
            print(f"  {action:8s} {size / 1e6:9.1f} MB  {key}  ->  {os.path.relpath(dest, a.root)}")
        if len(todo) > 20:
            print(f"  ... {len(todo) - 20} more")
        return 1 if a.check and any(n for v, (n, _) in verbs.items() if v not in CLEAN) else 0
    done, failed = 0, []
    try:
        with cf.ThreadPoolExecutor(max_workers=a.jobs) as ex:
            futs = {ex.submit(fetch, key, size, dest, endpoint, action): (key, size, etag, stamp)
                    for key, size, etag, stamp, dest, action in todo}
            for fut in cf.as_completed(futs):
                key, size, etag, stamp = futs[fut]
                try:
                    status = fut.result()
                except Exception as e:
                    failed.append((key, e))
                    print(f"  FAILED {key}: {e}", flush=True); continue
                if status != "OLDER":                 # nobody has vouched for those
                    manifest[product + key] = {"size": size, "etag": etag, "last_modified": stamp}
                done += size
                print(f"  {status:8s} {size / 1e6:9.1f} MB  {key}   [{done / 1e9:.2f} of {total / 1e9:.2f} GB]", flush=True)
    finally:
        # Written even when the pool is interrupted: what landed is recorded,
        # and what did not is judged by its age on the next run.
        save_manifest(a.root, manifest)
    if failed:
        # A caller that chains on success must not go on to submit a run
        # against a tree with holes in it, so the exit status carries this.
        print(f"{len(failed)} of {len(todo)} transfers FAILED:", flush=True)
        for key, e in failed:
            print(f"  {key}: {e}", flush=True)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
