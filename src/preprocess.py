#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
preprocess.py  — Calcinosis / CREST Paper 01
---------------------------------------------
QC and preprocessing for GSE195452 (Gur et al. 2022).

GSE195452 contains single-cell multiome (RNA + ATAC) data from SSc and
healthy skin, with 49 fibroblast subsets. For this paper we use only
the RNA modality.

Pipeline:
    1. Load all per-sample .h5 files from data/raw/GSE195452/
    2. Per-sample QC (mito %, gene count, doublet removal with Scrublet)
    3. Concatenate, normalise, HVG selection
    4. Save data/processed/calcinosis_qc.h5ad

Usage:
    python3 src/preprocess.py

Matches paper 01 (scleroderma-scvi) preprocess.py structure exactly.
Environment: ssc-scvi (see environment.yml)
"""

from pathlib import Path
import warnings
import gc

import numpy as np
import pandas as pd
import scanpy as sc
import anndata as an
from scipy import sparse
try:
    import scipy.sparse as sp
except Exception:
    from scipy import sparse as sp

warnings.filterwarnings("ignore", message=".*observed=False is deprecated.*")
warnings.filterwarnings("ignore", message=".*Variable names are not unique.*")
warnings.filterwarnings("ignore", message=".*Observation names are not unique.*")

RAW_DIR  = Path("data/raw/GSE138669")
PROC_DIR = Path("data/processed")
PER_DIR  = PROC_DIR / "per_sample"
PROC_DIR.mkdir(parents=True, exist_ok=True)
PER_DIR.mkdir(parents=True, exist_ok=True)

# QC thresholds — matched to paper 01 defaults
MIN_GENES   = 300
MAX_GENES   = 6000
MAX_PCT_MT  = 20.0
MIN_COUNTS  = 500
DOUBLET_THR = 0.25     # Scrublet score threshold


# ── helpers ──────────────────────────────────────────────────────────────────

def _strip_ens_version(idx):
    return pd.Index(idx).astype(str).str.replace(r"\.\d+$", "", regex=True)


def _symbol_series(var, fallback_index):
    for c in ("gene_symbol", "gene_symbols", "symbol", "SYMBOL",
              "gene_name", "features", "gene_names"):
        if c in var.columns:
            return var[c].astype(str)
    return pd.Index(fallback_index).astype(str)


def annotate_and_qc(ad):
    """Add gene symbols, mito/ribo flags, and compute QC metrics."""
    if "counts" not in ad.layers:
        ad.layers["counts"] = ad.X.copy()

    if "gene_symbol" not in ad.var.columns:
        ad.var["gene_symbol"] = _symbol_series(ad.var, ad.var_names)

    symu = ad.var["gene_symbol"].astype(str).str.upper()
    ad.var["mt"]   = symu.str.startswith(("MT-", "MT.", "MT_"))
    ad.var["ribo"] = symu.str.startswith(("RPS", "RPL"))
    sc.pp.calculate_qc_metrics(ad, qc_vars=["mt", "ribo"],
                               layer="counts", inplace=True)


def load_one_10x(h5_path: Path, sample_id: str) -> sc.AnnData:
    ad = sc.read_10x_h5(str(h5_path))

    if "gene_ids" in ad.var.columns:
        ens = _strip_ens_version(ad.var["gene_ids"])
    else:
        ens = _strip_ens_version(ad.var_names)
    ad.var_names = ens
    ad.var_names_make_unique()

    ad.obs_names = [f"{sample_id}_{bc}" for bc in ad.obs_names]
    ad.obs["sample"] = sample_id

    # Infer condition from sample name:
    # GSE195452 samples are named e.g. HC_1, SSc_2 etc.
    if False:  # GSE138669 is SSc-only
        ad.obs["condition"] = "Healthy"
    else:
        ad.obs["condition"] = "SSc"  # GSE138669: all samples are SSc

    annotate_and_qc(ad)
    return ad


def qc_filter(ad: sc.AnnData, sample_id: str) -> sc.AnnData:
    """Apply QC filters and optionally run Scrublet."""
    n0 = ad.n_obs

    # Basic numeric filters
    mask = (
        (ad.obs["n_genes_by_counts"] >= MIN_GENES) &
        (ad.obs["n_genes_by_counts"] <= MAX_GENES) &
        (ad.obs["total_counts"]      >= MIN_COUNTS) &
        (ad.obs["pct_counts_mt"]     <= MAX_PCT_MT)
    )
    ad = ad[mask].copy()
    print(f"  [{sample_id}] {n0} -> {ad.n_obs} cells after basic QC")

    # Scrublet doublet removal
    try:
        import scrublet as scr
        X = ad.layers["counts"]
        if sp.issparse(X):
            X = X.toarray()
        scrub = scr.Scrublet(X)
        scores, _ = scrub.scrub_doublets(verbose=False)
        ad.obs["doublet_score"] = scores
        doublet_mask = scores < DOUBLET_THR
        n_pre = ad.n_obs
        ad = ad[doublet_mask].copy()
        print(f"  [{sample_id}] {n_pre} -> {ad.n_obs} cells after Scrublet "
              f"(thr={DOUBLET_THR})")
    except ImportError:
        print(f"  [{sample_id}] Scrublet not installed — skipping doublet removal")
    except Exception as e:
        print(f"  [{sample_id}] Scrublet error ({e}) — skipping")

    return ad


def load_all_samples() -> sc.AnnData:
    """Load all .h5 files from RAW_DIR and concatenate."""
    h5_files = sorted(RAW_DIR.glob("*.h5"))
    if not h5_files:
        # Try subdirectories (common GEO layout)
        h5_files = sorted(RAW_DIR.rglob("*.h5"))
    if not h5_files:
        raise FileNotFoundError(
            f"No .h5 files found in {RAW_DIR}.\n"
            "Run: python3 src/data_download.py  then extract GSE195452_RAW.tar\n"
            "and place the per-sample .h5 files in data/raw/GSE195452/"
        )

    print(f"[info] Found {len(h5_files)} .h5 files in {RAW_DIR}")
    adatas = []
    for h5 in h5_files:
        # derive sample_id from filename, e.g. GSM5843975_HC_1_filtered_feature_bc_matrix.h5
        sample_id = h5.stem
        # strip common suffixes
        for suffix in ("raw_feature_bc_matrix", "_filtered_feature_bc_matrix", "_raw_feature_bc_matrix",
                       "_filtered", "_raw"):
            sample_id = sample_id.replace(suffix, "")
        # strip GSM prefix if present (keep condition label)
        if "_" in sample_id:
            parts = sample_id.split("_")
            # try to find HC/SSc label
            for i, p in enumerate(parts):
                if p.upper() in ("HC", "SSC", "HEALTHY", "DISEASE"):
                    sample_id = "_".join(parts[i:])
                    break

        print(f"  Loading {h5.name} as sample={sample_id}")
        try:
            ad = load_one_10x(h5, sample_id)
            ad = qc_filter(ad, sample_id)
            # cache per-sample h5ad
            out = PER_DIR / f"{sample_id}_qc.h5ad"
            ad.write_h5ad(out)
            adatas.append(ad)
            del ad
            gc.collect()
        except Exception as e:
            print(f"  [warn] Failed to load {h5.name}: {e}")

    if not adatas:
        raise RuntimeError("No samples loaded successfully.")

    print(f"\n[info] Concatenating {len(adatas)} samples...")
    adata = an.concat(adatas, join="outer", label="sample_batch",
                      keys=[a.obs["sample"].iloc[0] for a in adatas])
    adata.obs_names_make_unique()
    return adata


def preprocess(adata: sc.AnnData) -> sc.AnnData:
    """Normalise, log-transform, select HVGs."""
    print(f"[preprocess] {adata.n_obs} cells, {adata.n_vars} genes")

    # Ensure raw counts layer
    if "counts" not in adata.layers:
        adata.layers["counts"] = adata.X.copy()

    # Normalise to 10k counts per cell, log1p
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    # HVG selection — seurat_v3 on counts layer, corrected for batch
    sc.pp.highly_variable_genes(
        adata,
        flavor="seurat_v3",
        n_top_genes=4000,
        layer="counts",
        batch_key="sample" if "sample" in adata.obs.columns else None,
    )
    n_hvg = adata.var["highly_variable"].sum()
    print(f"[preprocess] {n_hvg} HVGs selected")

    # Scale (clip at 10 for stability)
    sc.pp.scale(adata, max_value=10)

    return adata


def main():
    print("=" * 60)
    print("Calcinosis / CREST Paper 01 — Preprocessing")
    print("Dataset: GSE195452 (Gur et al. 2022, Cell)")
    print("=" * 60)

    # Load and QC
    adata = load_all_samples()

    # Normalise + HVGs
    adata = preprocess(adata)

    # Save
    out = PROC_DIR / "calcinosis_qc.h5ad"
    adata.write_h5ad(out)
    print(f"\n[saved] {out}")
    print(f"  cells:  {adata.n_obs}")
    print(f"  genes:  {adata.n_vars}")
    print(f"  HVGs:   {adata.var['highly_variable'].sum()}")
    print(f"  obs:    {list(adata.obs.columns)}")


if __name__ == "__main__":
    main()
