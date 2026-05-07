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

DATA_DIR = Path(__file__).resolve().parent.parent / "data"
FORCING_DIR = DATA_DIR / "forcing"

GHUB_COLLECTION_ID = os.environ.get(
    "ISMIP7_GLOBUS_COLLECTION", "ccc9bbd2-4091-4e35-addd-eeb639cf5332"
)

OCEAN_BASE = "/ISMIP6/ISMIP7_Prep/AIS_ocean/share_with_modellers"

OCEAN_FILES = {
    "cesm2_waccm_historical": {
        "remote_dir": f"{OCEAN_BASE}/biascorr_cmip/CESM2-WACCM/historical/zhou_annual_30_sep/0yr/thetao_so_tf",
        "files": [
            "so_Oyr_CESM2-WACCM_historical_r1i1p1f1_ismip8km_60m_185001-201412.nc",
            "tf_Oyr_CESM2-WACCM_historical_r1i1p1f1_ismip8km_60m_185001-201412.nc",
            "thetao_Oyr_CESM2-WACCM_historical_r1i1p1f1_ismip8km_60m_185001-201412.nc",
        ],
        "local_dir": "ocean/cesm_waccm/historical",
    },
    "cesm2_waccm_ssp585": {
        "remote_dir": f"{OCEAN_BASE}/biascorr_cmip/CESM2-WACCM/ssp585/zhou_annual_30_sep/0yr/thetao_so_tf",
        "files": [
            "so_Oyr_CESM2-WACCM_ssp585_r1i1p1f1_ismip8km_60m_201501-210012.nc",
            "so_Oyr_CESM2-WACCM_ssp585_r1i1p1f1_ismip8km_60m_210101-215012.nc",
            "so_Oyr_CESM2-WACCM_ssp585_r1i1p1f1_ismip8km_60m_215101-220012.nc",
            "so_Oyr_CESM2-WACCM_ssp585_r1i1p1f1_ismip8km_60m_220101-225012.nc",
            "so_Oyr_CESM2-WACCM_ssp585_r1i1p1f1_ismip8km_60m_225101-229912.nc",
            "tf_Oyr_CESM2-WACCM_ssp585_r1i1p1f1_ismip8km_60m_201501-210012.nc",
            "tf_Oyr_CESM2-WACCM_ssp585_r1i1p1f1_ismip8km_60m_210101-215012.nc",
            "tf_Oyr_CESM2-WACCM_ssp585_r1i1p1f1_ismip8km_60m_215101-220012.nc",
            "tf_Oyr_CESM2-WACCM_ssp585_r1i1p1f1_ismip8km_60m_220101-225012.nc",
            "tf_Oyr_CESM2-WACCM_ssp585_r1i1p1f1_ismip8km_60m_225101-229912.nc",
            "thetao_Oyr_CESM2-WACCM_ssp585_r1i1p1f1_ismip8km_60m_201501-210012.nc",
            "thetao_Oyr_CESM2-WACCM_ssp585_r1i1p1f1_ismip8km_60m_210101-215012.nc",
            "thetao_Oyr_CESM2-WACCM_ssp585_r1i1p1f1_ismip8km_60m_215101-220012.nc",
            "thetao_Oyr_CESM2-WACCM_ssp585_r1i1p1f1_ismip8km_60m_220101-225012.nc",
            "thetao_Oyr_CESM2-WACCM_ssp585_r1i1p1f1_ismip8km_60m_225101-229912.nc",
        ],
        "local_dir": "ocean/cesm_waccm/ssp585",
    },
    "obs_climatology": {
        "remote_dir": f"{OCEAN_BASE}/climatology/zhou_annual_30_sep",
        "files": [
            "OI_Climatology_ismip8km_60m_so_extrap.nc",
            "OI_Climatology_ismip8km_60m_tf_extrap.nc",
            "OI_Climatology_ismip8km_60m_thetao_extrap.nc",
        ],
        "local_dir": "ocean/climatology",
    },
    "cesm2_waccm_climatology": {
        "remote_dir": f"{OCEAN_BASE}/biascorr_cmip/CESM2-WACCM/climatology",
        "files": [
            "ct_CESM2-WACCM_climatology_r1i1p1f1_ismip8km_60m_199501-201412.nc",
            "sa_CESM2-WACCM_climatology_r1i1p1f1_ismip8km_60m_199501-201412.nc",
        ],
        "local_dir": "ocean/cesm_waccm/climatology",
    },
    "cesm2_waccm_bias": {
        "remote_dir": f"{OCEAN_BASE}/biascorr_cmip/CESM2-WACCM/bias/zhou_annual_30_sep",
        "files": [
            "ct_CESM2-WACCM_bias_r1i1p1f1_ismip8km_60m_199501-201412.nc",
            "sa_CESM2-WACCM_bias_r1i1p1f1_ismip8km_60m_199501-201412.nc",
        ],
        "local_dir": "ocean/cesm_waccm/bias",
    },
}

CALIBRATION_FILES = {
    "observed_melt": {
        "remote_dir": f"{OCEAN_BASE}/meltmip",
        "files": [
            "melt_paolo_err_adusumilli_ismip8km.nc",
            "Melt_Paolo_Err_Adusumilli_imbie2.csv",
            "BFRN_ismip_8km.nc",
        ],
        "local_dir": "ocean/meltmip",
    },
    "ocean_model_data": {
        "remote_dir": f"{OCEAN_BASE}/meltmip/Ocean_Modelling_Data",
        "files": [
            "Mathiot23_cold_T.nc", "Mathiot23_cold_S.nc",
            "Mathiot23_cold_TF.nc", "Mathiot23_cold_m.nc",
            "Mathiod23_warm_T.nc", "Mathiod23_warm_S.nc",
            "Mathiod23_warm_TF.nc", "Mathiod23_warm_m.nc",
            "TimmermannUndGoeller2017_cold_T.nc", "TimmermannUndGoeller2017_cold_S.nc",
            "TimmermannUndGoeller2017_cold_TF.nc", "TimmermannUndGoeller2017_cold_m.nc",
            "TimmermannUndGoeller2017_warm_T.nc", "TimmermannUndGoeller2017_warm_S.nc",
            "TimmermannUndGoeller2017_warm_TF.nc", "TimmermannUndGoeller2017_warm_m.nc",
        ],
        "local_dir": "ocean/meltmip/ocean_modelling",
    },
    "basin_mask": {
        "remote_dir": f"{OCEAN_BASE}/imbie2",
        "files": ["basin_numbers_ismip8km.nc"],
        "local_dir": "ocean/imbie2",
    },
    "ismip_grid": {
        "remote_dir": f"{OCEAN_BASE}/ismip",
        "files": ["ismip_8km_60m_grid.nc"],
        "local_dir": "ocean/grid",
    },
    "topography": {
        "remote_dir": f"{OCEAN_BASE}/topography",
        "files": [
            "BedMachineAntarctica-v3_ismip_8km.nc",
            "BedMachineAntarctica-v3.nc",
            "Bedmap3_ismip_8km.nc",
            "bedmap3.nc",
        ],
        "local_dir": "ocean/topography",
    },
}


def get_globus_client():
    r"""Get an authenticated Globus TransferClient."""
    try:
        import globus_sdk
    except ImportError:
        print("ERROR: globus-sdk not installed.")
        print("  Install with: pip install globus-sdk")
        print("  Then run: python scripts/download_forcing.py --login")
        sys.exit(1)

    token_file = Path.home() / ".ismip7_globus_tokens.json"

    if token_file.exists():
        tokens = json.loads(token_file.read_text())
        authorizer = globus_sdk.AccessTokenAuthorizer(
            tokens.get("transfer_access_token", "")
        )
        client = globus_sdk.TransferClient(authorizer=authorizer)
        try:
            client.get_endpoint(GHUB_COLLECTION_ID)
            return client
        except globus_sdk.GlobusAPIError:
            print("  Stored token expired, re-authenticating...")

    return do_login()


def do_login():
    r"""Run the Globus native app auth flow."""
    import globus_sdk

    CLIENT_ID = "c9e8acfa-6c6d-4e68-aa6c-0ee4a12c4e2a"
    client = globus_sdk.NativeAppAuthClient(CLIENT_ID)
    client.oauth2_start_flow(
        requested_scopes="urn:globus:auth:scope:transfer.api.globus.org:all"
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
        "transfer_access_token": transfer_tokens["access_token"],
    }))
    token_file.chmod(0o600)
    print("  Tokens saved.")

    authorizer = globus_sdk.AccessTokenAuthorizer(
        transfer_tokens["access_token"]
    )
    return globus_sdk.TransferClient(authorizer=authorizer)


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
    print(f"\nDiscovering data at {OCEAN_BASE}/...")

    entries = list_remote_files(tc, OCEAN_BASE)
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

    remote_dir = file_set["remote_dir"]
    files = file_set["files"]
    local_dir = FORCING_DIR / file_set["local_dir"]
    local_dir.mkdir(parents=True, exist_ok=True)

    to_download = []
    for fn in files:
        local_path = local_dir / fn
        if local_path.exists():
            size_mb = local_path.stat().st_size / 1e6
            print(f"    {fn}  ({size_mb:.0f} MB) [exists]")
        else:
            to_download.append(fn)
            print(f"    {fn}  [needed]")

    if not to_download:
        return

    if dry_run:
        print(f"  Would download {len(to_download)} file(s)")
        return

    local_endpoint = os.environ.get("GLOBUS_LOCAL_ENDPOINT")
    if local_endpoint:
        td = globus_sdk.TransferData(
            tc, GHUB_COLLECTION_ID, local_endpoint,
            label=f"ISMIP7 {file_set['local_dir']}",
        )
        for fn in to_download:
            td.add_item(f"{remote_dir}/{fn}", str(local_dir / fn))
        result = tc.submit_transfer(td)
        task_id = result["task_id"]
        print(f"  Transfer submitted: {task_id}")
        print(f"  Monitor: https://app.globus.org/activity/{task_id}")
        return task_id

    print(f"\n  No GLOBUS_LOCAL_ENDPOINT set. To download these files:")
    print(f"  1. Install Globus Connect Personal: https://www.globus.org/globus-connect-personal")
    print(f"  2. Set GLOBUS_LOCAL_ENDPOINT=<your-endpoint-id>")
    print(f"  3. Re-run this script")
    print(f"  Or use the Globus web app to transfer from:")
    print(f"    {remote_dir}")
    print(f"  To your local machine at: {local_dir}")


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
        local_dir = FORCING_DIR / fset["local_dir"]
        total = len(fset["files"])
        if local_dir.exists():
            existing = [f for f in fset["files"] if (local_dir / f).exists()]
            if existing:
                total_mb = sum((local_dir / f).stat().st_size for f in existing) / 1e6
                print(f"  {name:30s}  {len(existing)}/{total} files ({total_mb:.0f} MB)")
            else:
                print(f"  {name:30s}  0/{total} files")
        else:
            print(f"  {name:30s}  not downloaded")


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
