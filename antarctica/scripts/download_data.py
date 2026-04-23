#!/usr/bin/env python3
"""
Download Antarctic datasets for ISMIP7 submission.

Required datasets:
    1. BedMachine Antarctica v4 (NSIDC-0756) - bed, thickness, surface, mask
    2. MEaSUREs Antarctic Ice Velocity v2 (NSIDC-0484) - vx, vy, errors
    3. RACMO2.4p1 Antarctic SMB (Zenodo) - surface mass balance

Requirements:
    pip install earthaccess

    You also need a free NASA Earthdata account:
    https://urs.earthdata.nasa.gov/users/new

Usage:
    python download_data.py
"""

import os
import sys
import urllib.request
from pathlib import Path

DATA_DIR = Path(__file__).resolve().parent.parent / "data"


def progress_hook(count, block_size, total_size):
    """Show download progress."""
    if total_size > 0:
        pct = min(100, count * block_size * 100 / total_size)
        mb = count * block_size / 1e6
        total_mb = total_size / 1e6
        sys.stdout.write(f"\r  {pct:.1f}% ({mb:.0f}/{total_mb:.0f} MB)")
    else:
        mb = count * block_size / 1e6
        sys.stdout.write(f"\r  {mb:.0f} MB downloaded")
    sys.stdout.flush()


def download_bedmachine():
    """Download BedMachine Antarctica v4 from NSIDC.

    Contains: bed elevation, ice thickness, surface elevation, ice mask,
    firn correction, geoid, source, error estimates.
    Resolution: 500 m, CRS: EPSG:3031
    """
    import earthaccess

    out_dir = DATA_DIR / "bedmachine"
    out_dir.mkdir(parents=True, exist_ok=True)

    existing = list(out_dir.glob("*BedMachine*.nc")) + list(
        out_dir.glob("*NSIDC-0756*.nc")
    )
    if existing:
        print(f"  Already downloaded: {existing[0].name}")
        return

    print("  Searching NSIDC for BedMachine Antarctica v4 (NSIDC-0756)...")
    results = earthaccess.search_data(short_name="NSIDC-0756", version="4", count=5)
    if not results:
        print("  ERROR: Could not find BedMachine Antarctica v4.")
        print("  Try manually from: https://nsidc.org/data/nsidc-0756")
        return

    print(f"  Found {len(results)} granule(s). Downloading...")
    earthaccess.download(results, str(out_dir))
    print("  BedMachine download complete!")


def download_velocity():
    """Download MEaSUREs InSAR-Based Antarctic Ice Velocity Map v2 from NSIDC.

    Contains: VX, VY (velocity components), STDX, STDY (std dev),
    ERRX, ERRY (formal errors), CNT (measurement count).
    Resolution: 450 m, CRS: EPSG:3031
    Dataset: NSIDC-0484 v2
    """
    import earthaccess

    out_dir = DATA_DIR / "velocity"
    out_dir.mkdir(parents=True, exist_ok=True)

    existing = list(out_dir.glob("*.nc"))
    if existing:
        print(f"  Already downloaded: {existing[0].name}")
        return

    print("  Searching NSIDC for MEaSUREs Antarctic Velocity v2 (NSIDC-0484)...")
    results = earthaccess.search_data(short_name="NSIDC-0484", version="2", count=5)
    if not results:
        print("  ERROR: Could not find MEaSUREs velocity v2.")
        print("  Try manually from: https://nsidc.org/data/nsidc-0484")
        return

    print(f"  Found {len(results)} granule(s). Downloading...")
    earthaccess.download(results, str(out_dir))
    print("  Velocity download complete!")


def download_racmo():
    """Download RACMO2.4p1 Antarctic SMB from Zenodo (public, no auth needed).

    Contains: monthly surface mass balance (glaciated) 1979-2023,
    plus grid/ice masks.
    Resolution: 11 km (ANT11 domain)
    DOI: 10.5281/zenodo.14217231
    """
    out_dir = DATA_DIR / "racmo"
    out_dir.mkdir(parents=True, exist_ok=True)

    base_url = "https://zenodo.org/records/14217232/files"
    files = [
        ("ANT11_masks.nc", "Grid/ice masks"),
        ("smbgl_monthlyS_ANT11_RACMO2.4p1_ERA5_197901_202312.nc", "SMB (glaciated)"),
    ]

    for fn, desc in files:
        out_path = out_dir / fn
        if out_path.exists():
            print(f"  {desc}: already downloaded")
            continue

        url = f"{base_url}/{fn}?download=1"
        print(f"  Downloading {desc} ({fn})...")
        try:
            urllib.request.urlretrieve(url, str(out_path), reporthook=progress_hook)
            print()  # newline after progress bar
        except Exception as e:
            print(f"\n  ERROR downloading {fn}: {e}")
            if out_path.exists():
                out_path.unlink()

    print("  RACMO download complete!")


def main():
    print("=" * 60)
    print("Antarctic Data Download for ISMIP7")
    print("=" * 60)

    DATA_DIR.mkdir(exist_ok=True)

    # 1. RACMO (no auth needed, do first)
    print("\n[1/3] RACMO2.4p1 SMB (Zenodo - public, no auth)")
    download_racmo()

    # 2-3. NSIDC products (need Earthdata auth)
    print("\nSetting up NASA Earthdata authentication...")
    print("  (You need a free account: https://urs.earthdata.nasa.gov/users/new)")
    try:
        import earthaccess

        earthaccess.login(persist=True)
    except ImportError:
        print("\nERROR: 'earthaccess' package not installed.")
        print("  Install with: pip install earthaccess")
        print("  Then re-run this script.")
        sys.exit(1)
    except Exception as e:
        print(f"\nERROR during Earthdata login: {e}")
        print("  Make sure you have valid Earthdata credentials.")
        sys.exit(1)

    print("\n[2/3] BedMachine Antarctica v4 (NSIDC-0756)")
    download_bedmachine()

    print("\n[3/3] MEaSUREs Antarctic Ice Velocity v2 (NSIDC-0484)")
    download_velocity()

    print("\n" + "=" * 60)
    print("All downloads complete!")
    print(f"Data directory: {DATA_DIR.resolve()}")
    print()
    print("Contents:")
    for p in sorted(DATA_DIR.rglob("*.nc")):
        size_mb = p.stat().st_size / 1e6
        print(f"  {p.relative_to(DATA_DIR)}  ({size_mb:.0f} MB)")
    print("=" * 60)


if __name__ == "__main__":
    main()
