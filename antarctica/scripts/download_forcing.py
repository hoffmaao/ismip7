#!/usr/bin/env python3
r"""Download ISMIP7 Antarctic forcing data from Globus.

Usage:
    python scripts/download_forcing.py --login
    python scripts/download_forcing.py --list
    python scripts/download_forcing.py --status
    python scripts/download_forcing.py --ocean
    python scripts/download_forcing.py --calibration
    python scripts/download_forcing.py
"""

import os, sys, json, argparse
from pathlib import Path

# Downloads land directly in the runtime tree that icepack2_tools/forcing.py
# reads (ISMIP7/AIS/...), mirroring the remote layout — no rename step.
FORCING_DIR = Path(__file__).resolve().parents[2] / "ISMIP7" / "AIS"

GHUB_COLLECTION_ID = os.environ.get(
    "ISMIP7_GLOBUS_COLLECTION", "ccc9bbd2-4091-4e35-addd-eeb639cf5332"
)

# The GHub collection is being ACTIVELY reorganized. State as of 2026-07-08,
# CONFIRMED UNCHANGED on a 2026-07-18 re-walk (full recursive listing, 32
# files total; headless CLI-token auth still works; --status: all 8 sets
# below fully mirrored). Still absent upstream: every scenario dir
# (ssp126/370/585 for BOTH ESMs), MRI-ESM2-0 entirely, CMIPraw, per-year
# atmosphere/ocean anything. Cores 1-6/8 remain unsourceable until the
# forcing group republishes - ask them (and where the ssp585 SDBN1-8000m
# atm + ocean we hold locally originally came from).
# (re-walked the whole tree — see scripts/preflight.py for the needs side):
#   /ISMIP6/ISMIP7_Prep/CMIP6_test_protocol/{AIS, GrIS, Tools, test}
# and under AIS/ only:
#   CESM2-WACCM/{bias, climatology, historical/ocean/extras}  (climatologies
#     + bias-correction ingredients ONLY — NO per-year forcing, no scenario dirs)
#   grid, obs (mipkit / OI-climatology / IMBIE / topography)
# Everything listed in OCEAN_FILES/CALIBRATION_FILES below is present and
# ALREADY MIRRORED locally, so this script is currently a no-op.
#
# NOT on the collection right now (was here days ago; removed mid-reorg):
#   * ../CMIPraw/  — the entire raw-CMIP tree is GONE.
#   * processed per-year atmosphere (SDBN1-8000m) and ocean (tf/so/thetao)
#     for ssp585 or any scenario — the CESM2-WACCM/ssp585 tree we hold in
#     ISMIP7/AIS/ (133 GB, atmosphere+ocean+fracture) is NO LONGER on the
#     share: treat the local copy as IRREPLACEABLE and back it up.
#   * MRI-ESM2-0 anything (AIS/MRI-ESM2-0 does not exist).
# => cores 1-6/8/10 cannot be sourced here; the processed forcing must be
#    (re)published by the forcing group, or regenerated from raw CMIP with a
#    processing pipeline not present here (Tools/ has only the GrIS pipeline).
AIS_BASE = "/ISMIP6/ISMIP7_Prep/CMIP6_test_protocol/AIS"

_LTM = "{var}_CESM2-WACCM_ltm_SDBN1_1960-1989.nc"
OCEAN_FILES = {
    "sdbn1_ltm_1960_1989": {
        # The exact climatology the acabf-anomaly files are referenced to.
        "remote_dir": f"{AIS_BASE}/CESM2-WACCM/climatology/SDBN1",
        "per_var_dirs": "v1",
        "files": [_LTM.format(var=v) for v in
                  ("acabf", "precip", "runoff", "snowfall", "snowmelt",
                   "tas", "ts")],
        "local_dir": "CESM2-WACCM/climatology/SDBN1",
    },
    "cesm2_waccm_ocean_climatology": {
        "remote_dir": f"{AIS_BASE}/CESM2-WACCM/climatology",
        "files": [
            "ct_CESM2-WACCM_climatology_r1i1p1f1_ismip8km_60m_1850_199501-201412.nc",
            "sa_CESM2-WACCM_climatology_r1i1p1f1_ismip8km_60m_1850_199501-201412.nc",
        ],
        "local_dir": "CESM2-WACCM/climatology",
    },
    "cesm2_waccm_ocean_bias": {
        "remote_dir": f"{AIS_BASE}/CESM2-WACCM/bias/zhou_annual_30_sep",
        "files": [
            "ct_CESM2-WACCM_bias_r1i1p1f1_ismip8km_60m_1850_199501-201412.nc",
            "sa_CESM2-WACCM_bias_r1i1p1f1_ismip8km_60m_1850_199501-201412.nc",
        ],
        "local_dir": "CESM2-WACCM/bias/zhou_annual_30_sep",
    },
    "oi_climatology_06nov": {
        "remote_dir": f"{AIS_BASE}/obs/ocean/climatology/zhou_annual_06_nov",
        "per_var_dirs": "v3",
        "files": [
            f"{v}_AIS_obs_ocean_climatology_zhou_annual_06_nov_v3_1972-2024.nc"
            for v in ("tf", "so", "thetao")
        ],
        "local_dir": "obs/ocean/climatology/zhou_annual_06_nov",
    },
}

CALIBRATION_FILES = {
    "obs_mipkit": {
        "remote_dir": f"{AIS_BASE}/obs/mipkit",
        "files": ["AntarcticaObsISMIP7-v1.1.nc"],
        "local_dir": "obs/mipkit",
    },
    "imbie_basins_v3": {
        "remote_dir": f"{AIS_BASE}/obs/ocean/IMBIE-basins/v3",
        "files": ["IMBIE-basins_AIS_obs_ocean_v3.nc"],
        "local_dir": "obs/ocean/IMBIE-basins/v3",
    },
    "ismip_grid": {
        "remote_dir": f"{AIS_BASE}/grid/ocean/ISMIP7/8km-60m/v3",
        "files": ["ISMIP7_8km-60m_AIS_grid_ocean_v3.nc"],
        "local_dir": "grid/ocean/ISMIP7/8km-60m/v3",
    },
    "topography": {
        "remote_dir": f"{AIS_BASE}/obs/ocean/topography",
        "per_file_dirs": {
            "BedMachineAntarctica-v3_AIS_obs_ocean_topography_v3.nc":
                "BedMachineAntarctica-v3/v3",
            "bedmap3_AIS_obs_ocean_topography_v3.nc": "bedmap3/v3",
        },
        "files": [
            "BedMachineAntarctica-v3_AIS_obs_ocean_topography_v3.nc",
            "bedmap3_AIS_obs_ocean_topography_v3.nc",
        ],
        "local_dir": "obs/ocean/topography",
    },
}


def _remote_local_pair(file_set, fn):
    r"""(remote_path, local_relpath) for one file, honoring per-var/v layout."""
    remote_dir = file_set["remote_dir"]
    sub = ""
    if "per_var_dirs" in file_set:
        var = fn.split("_")[0]
        sub = f"{var}/{file_set['per_var_dirs']}"
    elif "per_file_dirs" in file_set:
        sub = file_set["per_file_dirs"][fn]
    rd = f"{remote_dir}/{sub}".rstrip("/")
    ld = Path(file_set["local_dir"]) / sub
    return f"{rd}/{fn}", ld / fn


def _client_from_cli_storage():
    r"""TransferClient from the globus CLI's own token store
    (~/.globus/cli/storage.db): the CLI keeps a long-lived refresh token
    for transfer.api plus its per-user confidential-client credentials, so
    a machine that has ever run `globus login` works headlessly."""
    import sqlite3
    import globus_sdk

    db_path = Path.home() / ".globus" / "cli" / "storage.db"
    if not db_path.exists():
        return None
    try:
        db = sqlite3.connect(str(db_path))
        cc_data = json.loads(db.execute(
            "SELECT config_data_json FROM config_storage "
            "WHERE config_name='auth_client_data'"
        ).fetchone()[0])
        tok = json.loads(db.execute(
            "SELECT token_data_json FROM token_storage "
            "WHERE resource_server='transfer.api.globus.org'"
        ).fetchone()[0])
        cc = globus_sdk.ConfidentialAppAuthClient(
            cc_data["client_id"], cc_data["client_secret"]
        )
        authorizer = globus_sdk.RefreshTokenAuthorizer(
            tok["refresh_token"], cc
        )
        return globus_sdk.TransferClient(authorizer=authorizer)
    except Exception:
        return None


def get_globus_client():
    r"""Authenticated TransferClient: own refresh-token cache, then the
    globus CLI's token store, then an interactive login."""
    try:
        import globus_sdk
    except ImportError:
        print("ERROR: globus-sdk not installed.")
        print("  Install with: pip install globus-sdk")
        sys.exit(1)

    token_file = Path.home() / ".ismip7_globus_tokens.json"
    if token_file.exists():
        tokens = json.loads(token_file.read_text())
        if "refresh_token" in tokens:
            auth_client = globus_sdk.NativeAppAuthClient(tokens["client_id"])
            authorizer = globus_sdk.RefreshTokenAuthorizer(
                tokens["refresh_token"], auth_client
            )
            client = globus_sdk.TransferClient(authorizer=authorizer)
            try:
                client.get_endpoint(GHUB_COLLECTION_ID)
                return client
            except globus_sdk.GlobusAPIError:
                print("  Cached refresh token invalid, trying CLI store...")

    client = _client_from_cli_storage()
    if client is not None:
        try:
            client.get_endpoint(GHUB_COLLECTION_ID)
            print("  Authenticated via the globus CLI token store.")
            return client
        except globus_sdk.GlobusAPIError:
            pass

    return do_login()


def do_login():
    r"""Globus native-app auth flow, storing a REFRESH token so this is a
    one-time step per machine."""
    import globus_sdk

    CLIENT_ID = "c9e8acfa-6c6d-4e68-aa6c-0ee4a12c4e2a"
    client = globus_sdk.NativeAppAuthClient(CLIENT_ID)
    client.oauth2_start_flow(
        requested_scopes=[
            "urn:globus:auth:scope:transfer.api.globus.org:all",
        ],
        refresh_tokens=True,
    )

    authorize_url = client.oauth2_get_authorize_url()
    print(f"\nPlease visit this URL to authenticate:")
    print(f"  {authorize_url}")
    print()
    auth_code = input("Enter the authorization code: ").strip()

    token_response = client.oauth2_exchange_code_for_tokens(auth_code)
    transfer_tokens = token_response.by_resource_server["transfer.api.globus.org"]

    token_file = Path.home() / ".ismip7_globus_tokens.json"
    token_file.write_text(json.dumps({
        "client_id": CLIENT_ID,
        "refresh_token": transfer_tokens["refresh_token"],
        "transfer_access_token": transfer_tokens["access_token"],
    }))
    token_file.chmod(0o600)
    print("  Tokens saved (refresh token: no more logins needed).")

    authorizer = globus_sdk.RefreshTokenAuthorizer(
        transfer_tokens["refresh_token"], client,
        access_token=transfer_tokens["access_token"],
        expires_at=transfer_tokens["expires_at_seconds"],
    )
    return globus_sdk.TransferClient(authorizer=authorizer)


def get_local_endpoint():
    r"""Destination endpoint: GLOBUS_LOCAL_ENDPOINT env, else the Globus
    Connect Personal id registered on this machine."""
    ep = os.environ.get("GLOBUS_LOCAL_ENDPOINT")
    if ep:
        return ep
    gcp_id = Path.home() / ".globusonline" / "lta" / "client-id.txt"
    if gcp_id.exists():
        return gcp_id.read_text().strip()
    return None


def list_remote_files(tc, path, recursive=False):
    try:
        entries = tc.operation_ls(GHUB_COLLECTION_ID, path=path)
        results = []
        for entry in entries:
            full_path = f"{path}/{entry['name']}"
            if entry["type"] == "dir" and recursive:
                results.extend(list_remote_files(tc, full_path, recursive=True))
            else:
                results.append({
                    "name": entry["name"],
                    "path": full_path,
                    "type": entry["type"],
                    "size": entry.get("size", 0),
                })
        return results
    except Exception as e:
        print(f"  Error listing {path}: {e}")
        return []


def discover_layout(tc):
    print(f"\nDiscovering data at {AIS_BASE}/...")

    entries = list_remote_files(tc, AIS_BASE)
    if not entries:
        print("  No files found. Check your Globus auth and endpoint.")
        return

    for entry in entries:
        print(f"  {entry['type']:4s}  {entry['name']}")
        if entry["type"] == "dir":
            sub_entries = list_remote_files(tc, entry["path"])
            for sub in sub_entries[:8]:
                size_mb = sub.get("size", 0) / 1e6
                suf = f"  ({size_mb:.0f} MB)" if sub["type"] == "file" else ""
                print(f"         {sub['type']:4s}  {sub['name']}{suf}")
            if len(sub_entries) > 8:
                print(f"         ... and {len(sub_entries) - 8} more")


def download_file_set(tc, file_set, dry_run=False):
    import globus_sdk

    to_download = []
    for fn in file_set["files"]:
        remote, rel = _remote_local_pair(file_set, fn)
        local_path = FORCING_DIR / rel
        if local_path.exists():
            size_mb = local_path.stat().st_size / 1e6
            print(f"    {fn}  ({size_mb:.0f} MB) [exists]")
        else:
            to_download.append((remote, local_path))
            print(f"    {fn}  [needed]")

    if not to_download:
        return

    if dry_run:
        print(f"  Would download {len(to_download)} file(s)")
        return

    local_endpoint = get_local_endpoint()
    if local_endpoint:
        td = globus_sdk.TransferData(
            source_endpoint=GHUB_COLLECTION_ID,
            destination_endpoint=local_endpoint,
            label=f"ISMIP7 {file_set['local_dir']}",
            verify_checksum=True, sync_level="checksum",
        )
        for remote, local_path in to_download:
            td.add_item(remote, str(local_path))
        result = tc.submit_transfer(td)
        task_id = result["task_id"]
        print(f"  Transfer submitted: {task_id}")
        print(f"  Monitor: https://app.globus.org/activity/{task_id}")
        return task_id

    print(f"\n  No local Globus endpoint found. Either:")
    print(f"  1. Install/start Globus Connect Personal, or")
    print(f"  2. Set GLOBUS_LOCAL_ENDPOINT=<endpoint-id>, or")
    print(f"  3. Transfer by hand in the Globus web app from:")
    print(f"     {file_set['remote_dir']}")


def download_ocean(tc, dry_run=False):
    print("\n" + "=" * 60)
    print("Ocean Forcing (CESM2-WACCM, 8km x 60m grid)")
    print("=" * 60)

    for name, fset in OCEAN_FILES.items():
        print(f"\n  [{name}]")
        print(f"  Remote: {fset['remote_dir']}")
        download_file_set(tc, fset, dry_run=dry_run)


def download_calibration(tc, dry_run=False):
    print("\n" + "=" * 60)
    print("Melt Calibration Data")
    print("=" * 60)

    for name, fset in CALIBRATION_FILES.items():
        print(f"\n  [{name}]")
        print(f"  Remote: {fset['remote_dir']}")
        download_file_set(tc, fset, dry_run=dry_run)


def print_status():
    print("\n" + "=" * 60)
    print("ISMIP7 Forcing Data Status")
    print("=" * 60)

    all_sets = {**OCEAN_FILES, **CALIBRATION_FILES}
    for name, fset in all_sets.items():
        total = len(fset["files"])
        existing = []
        for fn in fset["files"]:
            _, rel = _remote_local_pair(fset, fn)
            if (FORCING_DIR / rel).exists():
                existing.append(FORCING_DIR / rel)
        if existing:
            total_mb = sum(p.stat().st_size for p in existing) / 1e6
            print(f"  {name:30s}  {len(existing)}/{total} files ({total_mb:.0f} MB)")
        else:
            print(f"  {name:30s}  0/{total} files")


def main():
    parser = argparse.ArgumentParser(
        description="Download ISMIP7 AIS forcing data from Globus",
    )
    parser.add_argument("--login", action="store_true", help="Authenticate with Globus")
    parser.add_argument("--list", action="store_true", help="Discover remote directory layout")
    parser.add_argument("--status", action="store_true", help="Show local download status")
    parser.add_argument("--ocean", action="store_true", help="Download ocean forcing only")
    parser.add_argument("--calibration", action="store_true", help="Download calibration data only")
    parser.add_argument("--dry-run", action="store_true", help="Show what would be downloaded")
    args = parser.parse_args()

    if args.status:
        print_status()
        return

    if args.login:
        do_login()
        print("Login successful.")
        return

    tc = get_globus_client()

    if args.list:
        discover_layout(tc)
        return

    dry = args.dry_run

    if args.ocean:
        download_ocean(tc, dry_run=dry)
    elif args.calibration:
        download_calibration(tc, dry_run=dry)
    else:
        download_ocean(tc, dry_run=dry)
        download_calibration(tc, dry_run=dry)

    print_status()


if __name__ == "__main__":
    main()
