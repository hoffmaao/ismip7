r"""A re-sync sees a forcing file the focus groups replaced in place.

The ISMIP7 share replaces files under the same name and the same version: the
AIS OCX ``dacabfdz`` came up spatially shifted and was corrected in place
(discussion #45), and the MRI-ESM2-0 ssp126 ``lake_properties`` was a byte
copy of the CESM2-WACCM file until it was swapped (#41). A downloader that
skips on size, and an audit that compares version strings, see neither. The
mirror route keeps each object's ETag in a manifest, the audit reads that
manifest and judges a version by what the reader would open, and the Globus
route picks dotted versions and the renamed MRI product and can be told to
let the checksum sync look at complete groups.

Nothing here touches the network: the S3 listing, the object bodies and the
Globus listing are all faked.
"""
import io
import os
import sys

import pytest

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(REPO, "antarctica", "scripts"))

import audit_forcing_versions as audit   # noqa: E402
import download_forcing as globus_route  # noqa: E402
import download_mirror as mirror         # noqa: E402

LISTING = b"""<?xml version="1.0"?>
<ListBucketResult xmlns="http://s3.amazonaws.com/doc/2006-03-01/">
  <Contents>
    <Key>ismip7-ais-forcing/data/OCX/ocean/main/tf_AIS_OCX_ocean_main_v1_1950-2025.nc</Key>
    <LastModified>2026-08-29T12:12:19.000Z</LastModified>
    <ETag>"c985def5a8891f091ecae6bde5b6109e-40"</ETag>
    <Size>2658589916</Size>
  </Contents>
</ListBucketResult>"""

THEN = 1_000_000.0                       # the mirror wrote its object
BEFORE, AFTER = THEN - 3600.0, THEN + 3600.0


def _file(path, body=b"abcd", mtime=AFTER):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "wb") as f:
        f.write(body)
    os.utime(path, (mtime, mtime))
    return str(path)


class _Body(io.BytesIO):
    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.close()


def test_the_listing_carries_the_etag_and_the_stamp(monkeypatch):
    monkeypatch.setattr(mirror.urllib.request, "urlopen", lambda req, timeout=None: _Body(LISTING))
    (key, size, etag, stamp), = mirror.list_keys("data/OCX/", mirror.MIRROR + "ismip7-ais-forcing/", "ismip7-ais-forcing/")
    assert key == "data/OCX/ocean/main/tf_AIS_OCX_ocean_main_v1_1950-2025.nc"
    assert size == 2658589916 and etag == "c985def5a8891f091ecae6bde5b6109e-40"
    assert stamp == mirror._epoch("2026-08-29T12:12:19.000Z") > 0


def test_a_recorded_file_is_skipped_until_its_etag_moves(tmp_path):
    dest = _file(tmp_path / "f.nc")
    assert mirror.plan(4, "aaa", THEN, dest, {"etag": "aaa"}) == "skip"
    assert mirror.plan(4, "bbb", THEN, dest, {"etag": "aaa"}) == "REPLACED"
    # once whole and recorded, now the wrong length: the object changed
    assert mirror.plan(9, "bbb", THEN, dest, {"etag": "aaa"}) == "REPLACED"


def test_a_file_the_manifest_never_saw_is_judged_by_its_age(tmp_path):
    fresh = _file(tmp_path / "fresh.nc", mtime=AFTER)
    stale = _file(tmp_path / "stale.nc", mtime=BEFORE)
    assert mirror.plan(4, "aaa", THEN, fresh, None) == "adopt"
    # older than the mirror's object proves nothing: most of a Globus tree
    # predates the mirror itself, so it is reported, not refetched
    assert mirror.plan(4, "aaa", THEN, stale, None) == "OLDER"
    assert mirror.fetch("k", 4, stale, "http://unused/", "OLDER") == "OLDER"


def test_only_a_prefix_of_the_current_object_is_resumed(tmp_path):
    assert mirror.plan(4, "aaa", THEN, str(tmp_path / "absent.nc"), None) == "fetch"
    assert mirror.plan(9, "aaa", THEN, _file(tmp_path / "new.nc", mtime=AFTER), None) == "resume"
    assert mirror.plan(9, "aaa", THEN, _file(tmp_path / "old.nc", mtime=BEFORE), None) == "fetch"


def test_a_replacement_lands_whole_or_not_at_all(tmp_path, monkeypatch):
    dest = _file(tmp_path / "f.nc", body=b"old!")
    monkeypatch.setattr(mirror.urllib.request, "urlopen", lambda req, timeout=None: _Body(b"ne"))
    with pytest.raises(IOError, match="stale copy stays"):
        mirror.fetch("k", 4, dest, "http://unused/", "REPLACED")
    assert open(dest, "rb").read() == b"old!" and not os.path.exists(dest + ".new")

    monkeypatch.setattr(mirror.urllib.request, "urlopen", lambda req, timeout=None: _Body(b"new!"))
    assert mirror.fetch("k", 4, dest, "http://unused/", "REPLACED") == "REPLACED"
    assert open(dest, "rb").read() == b"new!" and not os.path.exists(dest + ".new")


def test_the_manifest_round_trips(tmp_path):
    assert mirror.load_manifest(str(tmp_path)) == {}
    mirror.save_manifest(str(tmp_path), {"p/k": {"size": 4, "etag": "aaa", "last_modified": THEN}})
    assert mirror.load_manifest(str(tmp_path))["p/k"]["etag"] == "aaa"


# --- the audit ------------------------------------------------------------

def _key(esm, scenario, product, variable, name):
    return "/".join(x for x in ("data", esm, scenario, product, variable, name) if x)


def _entry(root, key, body=b"abcd", etag="aaa", mtime=AFTER, on_disk=None):
    r"""A mirror listing row, and (unless ``on_disk`` is False) its local copy."""
    if on_disk is not False:
        _file(on_disk or mirror.local_path(str(root), key), body, mtime)
    return (key, len(body), etag, THEN)


def _status(root, entries, manifest=None):
    listing = lambda prefix: [e for e in entries if e[0].startswith(prefix)]   # noqa: E731
    rows = audit.audit(str(root), audit.mirror_entries(["CESM2-WACCM", "MRI-ESM2-0"], None, listing), manifest or {})
    return {row[0][1:]: (row[4], row[3], row[5]) for row in rows}


def test_current_means_what_the_reader_would_open(tmp_path):
    esm, atm = "CESM2-WACCM", "SDBN1-8000m"
    ok = _key(esm, "ssp585", atm, "acabf", f"acabf_AIS_{esm}_ssp585_{atm}_v2_2015.nc")
    behind = _key(esm, "ssp126", atm, "acabf", f"acabf_AIS_{esm}_ssp126_{atm}_v3_2015.nc")
    pinned = _key(esm, "ssp370", atm, "acabf", f"acabf_AIS_{esm}_ssp370_{atm}_v3_2015.nc")
    entries = [_entry(tmp_path, ok),
               _entry(tmp_path, behind, on_disk=False),
               _entry(tmp_path, pinned)]
    # ssp126: only the superseded v2 is here. ssp370: v3 is here, and so is
    # the v2 the reader asks for first, which is what a run would open.
    for scenario in ("ssp126", "ssp370"):
        _file(tmp_path / esm / scenario / atm / "acabf" / "v2" / f"acabf_AIS_{esm}_{scenario}_{atm}_v2_2015.nc")
    got = _status(tmp_path, entries)
    assert got[("ssp585", atm, "acabf")][:2] == ("ok", "v2")
    assert got[("ssp126", atm, "acabf")][0] == "BEHIND"
    assert got[("ssp370", atm, "acabf")][:2] == ("PINNED", "v2")


def test_a_replaced_file_and_an_older_one_are_told_apart(tmp_path):
    esm = "CESM2-WACCM"
    tf = _key(esm, "ssp585", "ocean", "tf", f"tf_AIS_{esm}_ssp585_ocean_v3_2015-2024.nc")
    so = _key(esm, "ssp585", "ocean", "so", f"so_AIS_{esm}_ssp585_ocean_v3_2015-2024.nc")
    entries = [_entry(tmp_path, tf, etag="new"), _entry(tmp_path, so, mtime=BEFORE)]
    manifest = {audit.PRODUCT + tf: {"size": 4, "etag": "old", "last_modified": BEFORE}}
    got = _status(tmp_path, entries, manifest)
    assert got[("ssp585", "ocean", "tf")][0] == "REPLACED"
    assert got[("ssp585", "ocean", "so")] == ("ok", "v3", 1)      # one older file, unproven, not a failure


def test_the_renamed_mri_product_and_the_dotted_fracture_version(tmp_path):
    esm = "MRI-ESM2-0"
    atm = _key(esm, "ssp585", "GEMB-SDBN1-8000m", "acabf", f"acabf_AIS_{esm}_ssp585_GEMB-SDBN1-8000m_v1_2015.nc")
    old_name = tmp_path / esm / "ssp585" / "SDBN1-8000m" / "acabf" / "v1" / os.path.basename(atm)
    frac = _key("CESM2-WACCM", "ssp585", "fracture", "", "ice_shelf_collapse_mask_cesm2waccm_ssp585_ismip7_8km-v2.1.nc")
    docx = _key(esm, "ssp585", "", "", "Atmospheric forcing README.docx")
    entries = [_entry(tmp_path, atm, on_disk=old_name), _entry(tmp_path, frac), _entry(tmp_path, docx, on_disk=False)]
    _file(tmp_path / "CESM2-WACCM" / "ssp585" / "fracture" / "v2" / "ice_shelf_collapse_mask_cesm2waccm_ssp585_ismip7_8km-v2.nc")
    got = _status(tmp_path, entries)
    assert got[("ssp585", "GEMB-SDBN1-8000m", "acabf")][0] == "ok (as SDBN1-8000m)"
    assert got[("ssp585", "fracture", "")][:2] == ("ok", "v2.1")   # the reader takes the highest
    assert ("ssp585", "", "") not in got                          # a .docx is no forcing product


def test_the_ocx_tree_is_audited_in_its_own_layout(tmp_path):
    r"""Seen on a cluster tree in September 2026: the OCX ``dacabfdz`` still at
    v1, the spatially shifted file of discussion #45, with v2 on the mirror."""
    src = "RACMO2.3p2-ERA"
    grad = f"data/OCX/{src}/SDBN1-8000m/dacabfdz/dacabfdz_AIS_{src}_OCX_SDBN1-8000m_v2_1979.nc"
    smb = f"data/OCX/{src}/SDBN1-8000m/acabf/acabf_AIS_{src}_OCX_SDBN1-8000m_v1_1979.nc"
    tf = "data/OCX/ocean/main/tf_AIS_OCX_ocean_main_v1_1950-2025.nc"
    entries = [_entry(tmp_path, grad, on_disk=False), _entry(tmp_path, smb), _entry(tmp_path, tf)]
    _file(tmp_path / "OCX" / src / "SDBN1-8000m" / "dacabfdz" / "v1" / f"dacabfdz_AIS_{src}_OCX_SDBN1-8000m_v1_1979.nc")
    listing = lambda prefix: [e for e in entries if e[0].startswith(prefix)]   # noqa: E731
    rows = audit.audit(str(tmp_path), audit.mirror_entries(["OCX"], None, listing), {})
    got = {row[0][1:]: (row[4], row[3]) for row in rows}
    assert got[(src, "SDBN1-8000m", "dacabfdz")] == ("BEHIND", "v1")
    assert got[(src, "SDBN1-8000m", "acabf")] == ("ok", "v1")       # pinned v2 absent, so v1 is what is read
    assert got[("ocean", "main", "")] == ("ok", "v1")
    # and the file lands where the readers look for it
    assert mirror.local_path("R", tf) == os.path.join("R", "OCX", "ocean", "main", "v1", os.path.basename(tf))


def test_extra_and_extras_are_found_below_the_row(tmp_path):
    r"""Seen on Quartz in September 2026, with the whole mirror on disk: eight
    ``extra`` and ``extras`` rows reading MISSING with every one of their files
    present. The mirror nests them below the level a product/variable row
    reaches, at ``<product>/extra/climatology/<variable>/<version>/``, so a
    search that stops at the top of the row finds no version at all."""
    esm, atm = "CESM2-WACCM", "SDBN1-8000m"
    extra = (f"data/{esm}/historical/{atm}/extra/climatology/acabf/"
             f"acabf_AIS_{esm}_historical_{atm}_v2_1960-1989.nc")
    extras = (f"data/{esm}/historical/ocean/extras/bias/tf/"
              f"tf_AIS_{esm}_historical_ocean_v3_1995-2014.nc")
    got = _status(tmp_path, [_entry(tmp_path, extra), _entry(tmp_path, extras)])
    assert got[("historical", atm, "extra")][:2] == ("ok", None)
    assert got[("historical", "ocean", "extras")][:2] == ("ok", None)
    assert audit.local_versions(str(tmp_path), esm, "historical", atm, "extra") == ["v2"]
    # a row whose versions do sit at the top is untouched by the deeper search
    plain = _key(esm, "ssp585", atm, "acabf", f"acabf_AIS_{esm}_ssp585_{atm}_v2_2015.nc")
    _entry(tmp_path, plain)
    assert audit.local_versions(str(tmp_path), esm, "ssp585", atm, "acabf") == ["v2"]


def test_versions_below_a_row_sort_by_number(tmp_path):
    r"""``v10`` is above ``v2``, and a plain string sort puts it below."""
    esm, atm = "MRI-ESM2-0", "GEMB-SDBN1-8000m"
    for v in ("v2", "v10"):
        _file(tmp_path / esm / "ssp585" / atm / "extra" / "climatology" / "pr" / v
              / f"pr_AIS_{esm}_ssp585_{atm}_{v}_1960-1989.nc")
    assert audit.local_versions(str(tmp_path), esm, "ssp585", atm, "extra") == ["v2", "v10"]


# --- the dry run as the completeness gate -----------------------------------

def _dry_run(monkeypatch, capsys, root, listing, *flags):
    def no_network(*args, **kwargs):
        raise AssertionError("a dry run opened a connection")
    monkeypatch.setattr(mirror, "list_keys", lambda prefix, endpoint, product: listing)
    monkeypatch.setattr(mirror.urllib.request, "urlopen", no_network)
    monkeypatch.setattr(sys, "argv", ["download_mirror.py", "--root", str(root), *flags, "data/CESM2-WACCM/"])
    return mirror.main(), capsys.readouterr().out


def test_a_dry_run_counts_the_plan_by_verb_and_writes_nothing(tmp_path, monkeypatch, capsys):
    r"""The audit's rows read ``ok`` over a tree with 2,755 files absent
    (issue #41), and the dry run's 20-line preview cannot be counted over
    94,000 objects. The plan by verb, and the rows with work under them, are
    what say whether every file is there and unchanged."""
    esm = "CESM2-WACCM"
    tf = _key(esm, "ssp585", "ocean", "tf", f"tf_AIS_{esm}_ssp585_ocean_v3_2015.nc")
    so = _key(esm, "ssp585", "ocean", "so", f"so_AIS_{esm}_ssp585_ocean_v3_2015.nc")
    lake = _key(esm, "ssp585", "fracture", "", "lake_properties_cesm2waccm_ssp585_ismip7_8km-v2.1.nc")
    listing = [_entry(tmp_path, tf), _entry(tmp_path, so, body=b"abcdef", on_disk=False),
               _entry(tmp_path, lake, etag="bbb")]
    mirror.save_manifest(str(tmp_path), {"ismip7-ais-forcing/" + k: {"size": 4, "etag": "aaa", "last_modified": THEN}
                                         for k in (tf, lake)})
    manifest = tmp_path / mirror.MANIFEST_NAME
    before = sorted(os.path.join(d, f) for d, _, fs in os.walk(tmp_path) for f in fs), manifest.read_bytes()

    verbs, rows = mirror.summarize([(k, size, etag, stamp, mirror.local_path(str(tmp_path), k), verb)
                                    for (k, size, etag, stamp), verb in zip(listing, ("skip", "fetch", "REPLACED"))])
    assert verbs == {"skip": [1, 4], "adopt": [0, 0], "fetch": [1, 6], "resume": [0, 0],
                     "REPLACED": [1, 4], "OLDER": [0, 0]}
    # a flat fracture directory stops at the product, as the audit's row does
    assert rows == {("fetch", f"{esm}/ssp585/ocean/so"): [1, 6], ("REPLACED", f"{esm}/ssp585/fracture"): [1, 4]}

    status, out = _dry_run(monkeypatch, capsys, tmp_path, listing, "--dry-run")
    assert status == 0
    assert "by verb: 1 skip (0.00 GB), 0 adopt, 1 fetch (0.00 GB), 0 resume, 1 REPLACED (0.00 GB), 0 OLDER" in out
    assert [line.split() for line in out.splitlines() if line.endswith(("ocean/so", "ssp585/fracture"))] == [
        ["REPLACED", "1", "files", "0.00", "GB", f"{esm}/ssp585/fracture"],
        ["fetch", "1", "files", "0.00", "GB", f"{esm}/ssp585/ocean/so"]]
    # --check is the same dry run with the verdict in its exit status
    assert _dry_run(monkeypatch, capsys, tmp_path, listing, "--check")[0] == 1
    assert _dry_run(monkeypatch, capsys, tmp_path, listing[:1], "--check")[0] == 0
    # nothing fetched, nothing swapped in, and the manifest byte for byte as it was
    assert before == (sorted(os.path.join(d, f) for d, _, fs in os.walk(tmp_path) for f in fs), manifest.read_bytes())


# --- the Globus route -------------------------------------------------------

def _share(tree):
    return lambda tc, path, recursive=False: [{"name": n, "type": t, "path": f"{path}/{n}", "size": 0}
                                              for n, t in tree.get(path, [])]


def test_the_globus_route_picks_dotted_versions_and_the_renamed_product(monkeypatch):
    base = f"{globus_route.ISMIP7_BASE}/MRI-ESM2-0/ssp585"
    monkeypatch.setattr(globus_route, "list_remote_files", _share({
        f"{base}/fracture": [("v2", "dir"), ("v2.1", "dir"), ("v10", "dir"), ("notes", "dir"), ("v3", "file")],
        base: [("GEMB-SDBN1-8000m", "dir"), ("GEMB-SDBN1-2000m", "dir"), ("ocean", "dir")],
    }))
    assert globus_route._pick_version(None, f"{base}/fracture") == "v10"
    assert globus_route._atmosphere_dir(None, base) == "GEMB-SDBN1-8000m"
    assert globus_route._pick_version(None, f"{base}/absent") is None
