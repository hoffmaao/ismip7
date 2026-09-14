#!/usr/bin/env python3
r"""Compare the forcing versions on disk with the ISMIP7 Source Cooperative
mirror, the data-freeze copy of record (discussions #37 and #40, Sep 2026).

    python antarctica/scripts/audit_forcing_versions.py [--root ISMIP7/AIS]
        [--esm CESM2-WACCM --esm MRI-ESM2-0] [--scenario ssp585 ...]

For every <ESM>/<scenario>/<product>/<variable> the mirror publishes, print
the mirror's version directories next to the ones under --root, and flag the
rows where the local copy is missing or behind. The mirror keeps only the
current version of each product, so "behind" means "must re-sync before the
production runs" (and the README must cite the version used).

The mirror is anonymous S3 over HTTPS. Two things the endpoint insists on:
prefixes are relative to the product (``data/<ESM>/...``, no ``AIS/`` level,
while the keys it returns carry the product name), and a browser-like
User-Agent (Python's default is refused with 403). One flat, paginated
listing per ESM; the two core ESMs take well under a minute.
"""
import argparse
import os
import re
import sys
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET

ENDPOINT = "https://data.source.coop/ismip/ismip7-ais-forcing/"
PRODUCT = "ismip7-ais-forcing/"
HEADERS = {"User-Agent": "curl/8"}


def list_keys(prefix):
    r"""Every key under ``prefix`` (flat listing, paginated)."""
    keys, token = [], None
    while True:
        url = ENDPOINT + "?list-type=2&max-keys=1000&prefix=" + urllib.parse.quote(prefix)
        if token:
            url += "&continuation-token=" + urllib.parse.quote(token)
        with urllib.request.urlopen(urllib.request.Request(url, headers=HEADERS), timeout=60) as r:
            root = ET.fromstring(r.read())
        ns = {"s3": root.tag.split("}")[0].strip("{")}
        keys += [k.find("s3:Key", ns).text[len(PRODUCT):] for k in root.findall("s3:Contents", ns)]
        nxt = root.find("s3:NextContinuationToken", ns)
        if nxt is None:
            return keys
        token = nxt.text


VERSION = re.compile(r"[_-](v\d+(?:\.\d+)*)(?:_|\.nc$)")   # _v2_ in most names, -v2.1.nc in fracture names


def mirror_versions(esms, scenarios):
    r"""{(esm, scenario, product, variable): {versions}} from the mirror.

    The mirror keeps no version directories: only the current version of
    each product is published and the version lives in the filename
    (``acabf_AIS_CESM2-WACCM_ssp585_SDBN1-8000m_v2_2015.nc``, fracture
    ``..._v2.1.nc``), so it is read from there.
    """
    out = {}
    for esm in esms:
        for key in list_keys(f"data/{esm}/"):
            parts = key.split("/")            # data, esm, scenario, [product], [variable], file
            if len(parts) < 4:
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
            m = VERSION.search(parts[-1])
            out.setdefault((esm, sc, pr, var), set()).add(m.group(1) if m else "?")
    return out


def local_versions(root, esm, scenario, product, variable):
    d = os.path.join(root, esm, scenario, product, variable) if variable else os.path.join(root, esm, scenario, product)
    if not os.path.isdir(d):
        return []
    vers = sorted(x for x in os.listdir(d) if os.path.isdir(os.path.join(d, x)) and x.startswith("v"))
    if vers:
        return vers
    found = set()
    for f in os.listdir(d):
        m = VERSION.search(f)
        if m:
            found.add(m.group(1))
    return sorted(found)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--root", default=os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))), "ISMIP7", "AIS"))
    ap.add_argument("--esm", action="append", default=None)
    ap.add_argument("--scenario", action="append", default=None)
    a = ap.parse_args()
    esms = a.esm or ["CESM2-WACCM", "MRI-ESM2-0"]
    mv = mirror_versions(esms, set(a.scenario) if a.scenario else None)
    behind = missing = 0
    print(f"{'ESM':12s} {'scenario':11s} {'product':18s} {'variable':16s} {'mirror':12s} {'local':12s} status")
    for (esm, sc, pr, v), vers in sorted(mv.items()):
        vers = sorted(vers)
        loc = local_versions(a.root, esm, sc, pr, v)
        if not loc:
            status = "MISSING"; missing += 1
        elif set(vers) & set(loc):
            status = "ok"
        else:
            status = "BEHIND"; behind += 1
        print(f"{esm:12s} {sc:11s} {pr:18s} {v:16s} {' '.join(vers):12s} {' '.join(loc) or '-':12s} {status}")
    print(f"\n{len(mv)} mirror entries: {missing} missing locally, {behind} behind (local version no longer on the mirror)")
    return 1 if behind else 0


if __name__ == "__main__":
    sys.exit(main())
