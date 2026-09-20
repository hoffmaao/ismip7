#!/usr/bin/env python3
r"""Compare the forcing on disk with the ISMIP7 Source Cooperative mirror, the
data-freeze copy of record (discussions #37 and #40, Sep 2026).

    python antarctica/scripts/audit_forcing_versions.py [--root ISMIP7/AIS]
        [--esm CESM2-WACCM --esm MRI-ESM2-0 --esm OCX] [--scenario ssp585 ...]

The default is the two core ESMs and the OCX tree, which core 11 reads and
which is laid out the other way round (``OCX/<source>/<product>/<variable>``
and ``OCX/ocean/<scenario>``): it is the tree discussion #45's shifted
``dacabfdz`` lives in, so an audit that skips it cannot see that case.

For every <ESM>/<scenario>/<product>/<variable> the mirror publishes, print
the mirror's versions next to the ones under --root and the one a run would
open, and flag the rows that need attention before the production runs:

    MISSING    nothing local
    BEHIND     the mirror's version is not on disk
    PINNED     it is on disk, but the reader resolves another one (it asks for
               ``forcing.ATMOSPHERE_VERSION`` / ``OCEAN_VERSION`` first, and
               takes the highest fracture version)
    REPLACED   right version, but the mirror's object has changed since it was
               fetched: same name, new content (discussions #45 and #41)

The mirror keeps only the current version of each product, so BEHIND, PINNED
and REPLACED all mean "re-sync before the production runs", and they set the
exit status. The README must cite the versions a run used. A file the
download manifest has never seen and which is older than the mirror's object
is counted as ``older``: unproven either way, see ``download_mirror.py --older``.

The mirror is anonymous S3 over HTTPS. Two things the endpoint insists on:
prefixes are relative to the product (``data/<ESM>/...``, no ``AIS/`` level,
while the keys it returns carry the product name), and a browser-like
User-Agent (Python's default is refused with 403). One flat, paginated
listing per ESM; the two core ESMs take well under a minute.
"""
import argparse
import os
import sys

_SCRIPTS = os.path.dirname(os.path.abspath(__file__))
_PROJECT = os.path.dirname(os.path.dirname(_SCRIPTS))
sys.path.insert(0, _PROJECT)
sys.path.insert(0, _SCRIPTS)

from download_mirror import (                                   # noqa: E402
    DEFAULT_PRODUCT, MIRROR, VERSION, list_keys, load_manifest, local_path, plan,
)
from icepack2_tools.forcing import (                            # noqa: E402
    ATMOSPHERE_PRODUCTS, ATMOSPHERE_VERSION, OCEAN_VERSION,
    _resolve_version, _version_subdirs,
)

PRODUCT = DEFAULT_PRODUCT + "/"
ENDPOINT = MIRROR + PRODUCT


def mirror_entries(esms, scenarios, listing=None):
    r"""{(esm, scenario, product, variable): [(key, size, etag, stamp), ...]}.

    The mirror keeps no version directories: only the current version of
    each product is published and the version lives in the filename
    (``acabf_AIS_CESM2-WACCM_ssp585_SDBN1-8000m_v2_2015.nc``, fracture
    ``..._v2.1.nc``), so it is read from there. Anything that is not NetCDF
    (the ``Atmospheric forcing README.docx`` of discussion #41) is no forcing
    product and carries no version, so it is no row.
    """
    listing = listing or (lambda prefix: list_keys(prefix, ENDPOINT, PRODUCT))
    out = {}
    for esm in esms:
        for entry in listing(f"data/{esm}/"):
            parts = entry[0].split("/")       # data, esm, scenario, [product], [variable], file
            if len(parts) < 4 or not parts[-1].endswith(".nc"):
                continue
            sc = parts[2]
            if scenarios and sc not in scenarios:
                continue
            if len(parts) == 4:               # data/esm/scenario/file (loose file)
                pr, var = "", ""
            elif len(parts) == 5:             # data/esm/scenario/fracture/file
                pr, var = parts[3], ""
            else:                             # data/esm/scenario/product/variable/file
                pr, var = parts[3], parts[4]
            out.setdefault((esm, sc, pr, var), []).append(entry)
    return out


def versions_of(entries):
    found = set()
    for key, *_ in entries:
        m = VERSION.search(os.path.basename(key))
        found.add(m.group(1) if m else "?")
    return sorted(found)


def local_product(root, esm, scenario, product):
    r"""The directory on disk that holds ``product``. MRI-ESM2-0's ``SDBN1-*``
    became ``GEMB-SDBN1-*`` in August 2026 with the data unchanged (#37), and
    the reader takes whichever is there, so a tree fetched before the rename
    is current under its old name."""
    if os.path.isdir(os.path.join(root, esm, scenario, product)):
        return product
    for name in ATMOSPHERE_PRODUCTS:
        if product.startswith(name + "-"):
            for other in ATMOSPHERE_PRODUCTS:
                alias = other + product[len(name):]
                if os.path.isdir(os.path.join(root, esm, scenario, alias)):
                    return alias
    return product


def local_versions(root, esm, scenario, product, variable):
    d = os.path.join(root, esm, scenario, product, variable) if variable else os.path.join(root, esm, scenario, product)
    if not os.path.isdir(d):
        return []
    vers = [name for _, name in _version_subdirs(d)]
    if vers:
        return vers
    found = set()
    for f in os.listdir(d):
        m = VERSION.search(f)
        if m:
            found.add(m.group(1))
    return sorted(found)


def resolved_version(root, esm, scenario, product, variable):
    r"""The version a run opens, by the readers' own rules, or None where no
    reader pins one (loose files, a flat fracture directory)."""
    d = os.path.join(root, esm, scenario, product, variable) if variable else os.path.join(root, esm, scenario, product)
    subdirs = _version_subdirs(d)
    if not subdirs:
        return None
    if product == "ocean":
        return _resolve_version(d, OCEAN_VERSION)
    if any(product.startswith(name + "-") for name in ATMOSPHERE_PRODUCTS):
        return _resolve_version(d, ATMOSPHERE_VERSION)
    return subdirs[-1][1]                     # fracture, and anything unpinned: the highest


def audit(root, entries, manifest):
    r"""``[(row, mirror versions, local versions, resolved, status, n_older)]``."""
    rows = []
    for (esm, sc, pr, v), keys in sorted(entries.items()):
        vers = versions_of(keys)
        on_disk = local_product(root, esm, sc, pr)
        loc = local_versions(root, esm, sc, on_disk, v)
        resolved = resolved_version(root, esm, sc, on_disk, v)
        verdicts = []
        for key, size, etag, stamp in keys:
            dest = local_path(root, key)
            if on_disk != pr:
                dest = dest.replace(os.sep + pr + os.sep, os.sep + on_disk + os.sep, 1)
            verdicts.append(plan(size, etag, stamp, dest, manifest.get(PRODUCT + key)))
        if not loc:
            status = "MISSING"
        elif not set(vers) & set(loc):
            status = "BEHIND"
        elif resolved is not None and resolved not in vers:
            status = "PINNED"
        elif "REPLACED" in verdicts:
            status = "REPLACED"
        else:
            status = "ok"
        if status == "ok" and on_disk != pr:
            status = f"ok (as {on_disk})"
        rows.append(((esm, sc, pr, v), vers, loc, resolved, status, verdicts.count("OLDER")))
    return rows


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--root", default=os.path.join(_PROJECT, "ISMIP7", "AIS"))
    ap.add_argument("--esm", action="append", default=None)
    ap.add_argument("--scenario", action="append", default=None)
    a = ap.parse_args()
    esms = a.esm or ["CESM2-WACCM", "MRI-ESM2-0", "OCX"]
    entries = mirror_entries(esms, set(a.scenario) if a.scenario else None)
    rows = audit(a.root, entries, load_manifest(a.root))
    print(f"{'ESM':12s} {'scenario':11s} {'product':18s} {'variable':16s} {'mirror':10s} {'local':10s} {'reads':6s} status")
    for (esm, sc, pr, v), vers, loc, resolved, status, older in rows:
        note = f"  ({older} older than the mirror's object)" if older else ""
        print(f"{esm:12s} {sc:11s} {pr:18s} {v:16s} {' '.join(vers):10s} {' '.join(loc) or '-':10s} "
              f"{resolved or '-':6s} {status}{note}")
    count = {s: sum(1 for r in rows if r[4] == s) for s in ("MISSING", "BEHIND", "PINNED", "REPLACED")}
    older = sum(r[5] for r in rows)
    print(f"\n{len(rows)} mirror entries: {count['MISSING']} missing locally, "
          f"{count['BEHIND']} behind (local version no longer on the mirror), "
          f"{count['PINNED']} pinned (the reader opens a version the mirror dropped), "
          f"{count['REPLACED']} replaced (same name, new content)")
    if older:
        print(f"{older} file(s) predate the download manifest and are older than the mirror's "
              f"object: unproven, see download_mirror.py --older")
    return 1 if count["BEHIND"] or count["PINNED"] or count["REPLACED"] else 0


if __name__ == "__main__":
    sys.exit(main())
