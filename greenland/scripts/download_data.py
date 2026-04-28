#! /usr/bin/env python3
# vim:fenc=utf-8
#
# Copyright © 2026 David Lilien <dlilien@iu.edu>
#
"""
Download Greenland datasets for ISMIP7 submission.

Required datasets:
    1. BedMachine Greenland v6 (IDBMG4) - bed, thickness, surface, mask
    2. MEaSUREs Greenland Ice Velocity
       - Multi-year Ice Sheet Velocity Mosaic v1 (NSIDC-0670) - vx, vy
       - Optional Annual Ice Velocity Mosaics v5 (NSIDC-0725) - vx, vy
    3. RACMO2.4p1 Greenland SMB (Zenodo) - surface mass balance
    4. PROMICE-2022 Ice Mask (GEUS Dataverse) - ice-sheet outline

Requirements:
    pip install earthaccess

    You also need a free NASA Earthdata account:
    https://urs.earthdata.nasa.gov/users/new

Usage:
    python download_data.py
"""

import argparse
import sys
import urllib.request
import zipfile
from pathlib import Path

DEFAULT_DATA_DIR = Path(__file__).resolve().parent.parent / "data"
PROMICE_ICE_MASK_FILENAME = "02-PROMICE-2022-IceMask-polygon-v3.gpkg"
PROMICE_ICE_MASK_URL = (
    "https://dataverse.geus.dk/api/access/datafile/:persistentId"
    "?persistentId=doi:10.22008/FK2/O8CLRE/NSEY7V&version=3.2"
)


def parse_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Download Greenland datasets for ISMIP7 submission."
    )
    parser.add_argument(
        "--data-dir",
        type=Path,
        default=DEFAULT_DATA_DIR,
        help="Directory for downloaded data. Defaults to data in the parent folder.",
    )
    parser.add_argument(
        "--include-annual-velocity",
        action="store_true",
        help="Also download MEaSUREs Greenland annual velocity mosaics.",
    )
    return parser.parse_args()


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


def download_bedmachine(data_dir):
    """Download IceBridge BedMachine Greenland v6 from NSIDC.

    Contains: bed elevation, ice thickness, surface elevation, ice/ocean/land
    mask, source and error estimates.
    Resolution: 150 m, CRS: EPSG:3413
    Dataset: IDBMG4 v6
    """
    import earthaccess

    out_dir = data_dir / "bedmachine"
    out_dir.mkdir(parents=True, exist_ok=True)

    existing = list(out_dir.glob("*BedMachineGreenland*.nc")) + list(
        out_dir.glob("*IDBMG4*.nc")
    )
    if existing:
        print(f"  Already downloaded: {existing[0].name}")
        return

    print("  Searching NSIDC for BedMachine Greenland v6 (IDBMG4)...")
    results = earthaccess.search_data(short_name="IDBMG4", version="6", count=10)
    if not results:
        print("  ERROR: Could not find BedMachine Greenland v6.")
        print("  Try manually from: https://nsidc.org/data/idbmg4/versions/6")
        return

    print(f"  Found {len(results)} granule(s). Downloading...")
    earthaccess.download(results, str(out_dir))
    print("  BedMachine download complete!")


def download_velocity(data_dir, include_annual=False):
    """Download MEaSUREs Greenland annual and multiyear velocity products.

    Multiyear product: NSIDC-0670 v1, 250 m, GeoTIFF, EPSG:3413
    Optional annual product: NSIDC-0725 v5, 200 m, GeoTIFF, EPSG:3413
    """
    import earthaccess

    products = [
        {
            "label": "MEaSUREs Greenland multiyear velocity v1",
            "short_name": "NSIDC-0670",
            "version": "1",
            "count": 10,
            "out_dir": data_dir / "velocity" / "multiyear",
            "manual_url": "https://nsidc.org/data/nsidc-0670/versions/1",
        },
    ]
    if include_annual:
        products.append(
            {
                "label": "MEaSUREs Greenland annual velocity v5",
                "short_name": "NSIDC-0725",
                "version": "5",
                "count": 100,
                "out_dir": data_dir / "velocity" / "annual",
                "manual_url": "https://nsidc.org/data/nsidc-0725/versions/5",
            }
        )

    for product in products:
        out_dir = product["out_dir"]
        out_dir.mkdir(parents=True, exist_ok=True)

        existing = list(out_dir.glob("*.tif")) + list(out_dir.glob("*.tiff"))
        if existing:
            print(f"  {product['label']}: already downloaded ({existing[0].name})")
            continue

        print(f"  Searching NSIDC for {product['label']} ({product['short_name']})...")
        results = earthaccess.search_data(
            short_name=product["short_name"],
            version=product["version"],
            count=product["count"],
        )
        if not results:
            print(f"  ERROR: Could not find {product['label']}.")
            print(f"  Try manually from: {product['manual_url']}")
            continue

        print(f"  Found {len(results)} granule(s). Downloading...")
        earthaccess.download(results, str(out_dir))
    print("  Velocity download complete!")


def download_racmo(data_dir):
    """Download RACMO2.4p1 Greenland data from Zenodo.

    Contains: monthly Greenland 11 km RACMO2.4p1 fields for 2006-2015,
    including surface mass balance, melt, refreezing, precipitation, runoff,
    SEB components, near-surface temperature, and wind speed.
    Resolution: 11 km (GRN11 domain)
    DOI: 10.5281/zenodo.13773130
    """
    out_dir = data_dir / "racmo"
    out_dir.mkdir(parents=True, exist_ok=True)

    fn = "GRN11_RACMO24p1_data.zip"
    zip_path = out_dir / fn
    extracted_marker = out_dir / ".GRN11_RACMO24p1_data_extracted"

    if extracted_marker.exists():
        print("  RACMO Greenland data: already downloaded and extracted")
        return

    if not zip_path.exists():
        url = f"https://zenodo.org/records/13773130/files/{fn}?download=1"
        print(f"  Downloading RACMO2.4p1 Greenland data ({fn})...")
        try:
            urllib.request.urlretrieve(url, str(zip_path), reporthook=progress_hook)
            print()  # newline after progress bar
        except Exception as e:
            print(f"\n  ERROR downloading {fn}: {e}")
            if zip_path.exists():
                zip_path.unlink()
            return
    else:
        print(f"  {fn}: already downloaded")

    print("  Extracting RACMO Greenland data...")
    with zipfile.ZipFile(zip_path) as zf:
        zf.extractall(out_dir)
    extracted_marker.touch()
    print("  RACMO download complete!")


def download_promice_ice_mask(data_dir):
    """Download the PROMICE-2022 Greenland Ice Mask from GEUS Dataverse.

    Contains: Greenland Ice Sheet outline from August 2022 as polygon vector
    features.
    Format: GeoPackage
    DOI: 10.22008/FK2/O8CLRE/NSEY7V
    """
    out_dir = data_dir / "promice"
    out_dir.mkdir(parents=True, exist_ok=True)

    out_path = out_dir / PROMICE_ICE_MASK_FILENAME
    if out_path.exists():
        print(f"  PROMICE ice mask: already downloaded ({out_path.name})")
        return

    print(f"  Downloading PROMICE ice mask ({PROMICE_ICE_MASK_FILENAME})...")
    try:
        urllib.request.urlretrieve(
            PROMICE_ICE_MASK_URL, str(out_path), reporthook=progress_hook
        )
        print()  # newline after progress bar
    except Exception as e:
        print(f"\n  ERROR downloading {PROMICE_ICE_MASK_FILENAME}: {e}")
        if out_path.exists():
            out_path.unlink()
        return

    print("  PROMICE ice mask download complete!")


def main():
    """Download all configured Greenland datasets from command-line arguments."""
    args = parse_args()
    data_dir = args.data_dir.expanduser()

    print("=" * 60)
    print("Greenland Data Download for ISMIP7")
    print("=" * 60)

    data_dir.mkdir(parents=True, exist_ok=True)

    # 1-2. Public products (no auth needed, do first)
    print("\n[1/4] RACMO2.4p1 Greenland SMB (Zenodo - public, no auth)")
    download_racmo(data_dir)

    print("\n[2/4] PROMICE-2022 Ice Mask (GEUS Dataverse - public, no auth)")
    download_promice_ice_mask(data_dir)

    # 3-4. NSIDC products (need Earthdata auth)
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

    print("\n[3/4] BedMachine Greenland v6 (IDBMG4)")
    download_bedmachine(data_dir)

    velocity_label = (
        "multiyear + annual" if args.include_annual_velocity else "multiyear"
    )
    print(f"\n[4/4] MEaSUREs Greenland Ice Velocity ({velocity_label})")
    download_velocity(data_dir, include_annual=args.include_annual_velocity)

    print("\n" + "=" * 60)
    print("All downloads complete!")
    print(f"Data directory: {data_dir.resolve()}")
    print()
    print("Contents:")
    for p in sorted(data_dir.rglob("*")):
        if not p.is_file() or p.name.startswith("."):
            continue
        size_mb = p.stat().st_size / 1e6
        print(f"  {p.relative_to(data_dir)}  ({size_mb:.0f} MB)")
    print("=" * 60)


if __name__ == "__main__":
    main()
