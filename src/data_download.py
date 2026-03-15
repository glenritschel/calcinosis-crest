#!/usr/bin/env python3
"""
data_download.py
----------------
Download GSE195452 (Gur et al. 2022, Cell) from GEO.

GSE195452 is an LGR5+ fibroblast atlas of SSc skin with high-resolution
fibroblast subsets. It is the primary dataset for the calcinosis paper
because calcinosis is driven by osteogenic/mineralising fibroblast states,
and this atlas resolves those states better than GSE138669.

Usage:
    python3 src/data_download.py

Outputs (into data/raw/GSE195452/):
    *.h5 or *_matrix.mtx.gz  (raw 10x counts)
    metadata files

Note: After download, copy the .h5 files into data/raw/GSE195452/
and run preprocess.py
"""

import argparse
from pathlib import Path
import urllib.request
import GEOparse

RAW = Path("data/raw")

GEO_ACCESSION = "GSE138669"

# Direct supplementary file URLs (fallback if GEOparse WAF blocks)
DIRECT_URLS = [
    (
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE195nnn/GSE195452/suppl/"
        "GSE195452_RAW.tar",
        "GSE195452_RAW.tar",
    ),
]


def download_geo(geo_accession: str = GEO_ACCESSION):
    dest = RAW / geo_accession
    dest.mkdir(parents=True, exist_ok=True)
    print(f"[info] downloading {geo_accession} into {dest.resolve()}")
    try:
        gse = GEOparse.get_GEO(geo=geo_accession, destdir=str(RAW), silent=False)
        print("[done] GEOparse download complete.")
        print("[info] Check data/raw/GSE195452 for supplementary .h5 / mtx files.")
        print(f"[info] Found {len(gse.gsms)} samples.")
    except Exception as exc:
        print(f"[warn] GEOparse failed: {exc}")
        print("[info] Trying direct FTP download...")
        _direct_download(dest)


def _direct_download(dest: Path):
    """Fallback: direct urllib download from NCBI FTP."""
    headers = {"User-Agent": "Mozilla/5.0 (compatible; research downloader)"}
    for url, fname in DIRECT_URLS:
        out = dest / fname
        if out.exists():
            print(f"[skip] {fname} already present.")
            continue
        print(f"[download] {url}")
        req = urllib.request.Request(url, headers=headers)
        with urllib.request.urlopen(req) as resp, open(out, "wb") as fh:
            while chunk := resp.read(1 << 20):
                fh.write(chunk)
        print(f"[done] {out}")
        if fname.endswith(".tar"):
            import tarfile
            print(f"[extract] {out}")
            with tarfile.open(out) as tf:
                tf.extractall(dest)
            print("[done] extraction complete.")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--geo", type=str, default=GEO_ACCESSION,
                    help="GEO accession (default: GSE195452)")
    args = ap.parse_args()
    download_geo(args.geo)
