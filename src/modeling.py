#!/usr/bin/env python3
"""
modeling.py  — Calcinosis / CREST Paper 01
-------------------------------------------
scVI embedding, UMAP, and Leiden clustering for GSE195452.

Pipeline:
    1. Load calcinosis_qc.h5ad
    2. Train scVI (batch_key="sample") on HVGs
    3. Compute UMAP and Leiden clusters
    4. Coarse cell-type annotation using canonical marker genes
    5. Save calcinosis_scvi.h5ad  and  calcinosis_scvi_annot.h5ad

Usage:
    python3 src/modeling.py

Matches paper 01 (scleroderma-scvi) modeling.py structure exactly.
Environment: ssc-scvi  (Python 3.10, scanpy 1.11.4, scvi-tools 1.3.3)
"""

from pathlib import Path
import numpy as np
import pandas as pd
import scanpy as sc
import scvi
import scipy.sparse as sp
import argparse

PROC  = Path("data/processed")
FIGS  = Path("figures");  FIGS.mkdir(parents=True, exist_ok=True)

# ── scVI hyperparameters (matched to paper 01) ───────────────────────────────
N_TOP_HVGS  = 4000
N_LATENT    = 30
N_LAYERS    = 2
N_HIDDEN    = 128
MAX_EPOCHS  = 400
BATCH_SIZE  = 256
EARLY_STOP  = True
LEIDEN_RES  = 0.5   # coarser than paper 01 to account for smaller dataset

SEEDS = {"numpy": 0, "torch": 0, "scvi": 0}

# ── Canonical marker genes for coarse annotation ─────────────────────────────
# Based on GSE195452 original cell-type labels (Gur et al. 2022)
# and calcinosis-relevant biology
MARKERS = {
    "Fibroblasts":          ["DCN", "LUM", "COL1A1", "COL3A1", "PDGFRA"],
    "LGR5+ Fibroblasts":    ["LGR5", "COMP", "MYOC", "IGFBP2"],
    "Myofibroblasts":       ["ACTA2", "TAGLN", "MYH11", "SFRP2", "PRSS23"],
    "Pericyte/SMC":         ["RGS5", "MCAM", "PDGFRB", "NOTCH3"],
    "Endothelial":          ["PECAM1", "CDH5", "VWF", "CLDN5"],
    "Lymphatic endothelium":["LYVE1", "PROX1", "PDPN"],
    "Keratinocytes (basal)":["KRT5", "KRT14", "TP63", "COL17A1"],
    "Keratinocytes (suprabasal)": ["KRT1", "KRT10", "FLG", "LOR"],
    "T cells":              ["CD3D", "CD3E", "CD3G", "TRAC"],
    "B cells":              ["CD19", "MS4A1", "CD79A"],
    "Plasma cells":         ["MZB1", "IGHG1", "JCHAIN", "XBP1"],
    "Myeloid (mono/mac)":   ["CD68", "CD14", "CSF1R", "MRC1"],
    "Mast cells":           ["TPSAB1", "CPA3", "HPGDS"],
    "NK cells":             ["NCAM1", "GNLY", "NKG7"],
}

# Calcinosis-specific gene signatures for scoring
CALCINOSIS_SIGNATURES = {
    "Osteogenic":     ["RUNX2", "SP7", "BGLAP", "ALPL", "COL10A1", "IBSP"],
    "Hypoxia":        ["HIF1A", "VEGFA", "LDHA", "ENO1", "SLC2A1", "BNIP3"],
    "Calcification":  ["ENPP1", "ANKH", "MGP", "FETUB", "TNAP", "ALPL",
                       "BMP2", "BMP4", "SMAD1", "SMAD5"],
    "Fibrosis":       ["TGFB1", "TGFB2", "CTGF", "FN1", "POSTN", "THBS1"],
    "Inflammation":   ["IL6", "IL1B", "TNF", "CXCL8", "CXCL12", "CCL2"],
}


def _symbol_series(a: sc.AnnData) -> pd.Series:
    for c in ["gene_symbol", "gene_symbols", "symbol", "SYMBOL",
              "gene_name", "gene_names", "feature_name", "features"]:
        if c in a.var.columns:
            return a.var[c].astype(str)
    return a.var_names.astype(str)


def _ensure_hvgs(adata, n_top_genes=N_TOP_HVGS, batch_key=None):
    if ("highly_variable" in adata.var.columns and
            adata.var["highly_variable"].dtype == bool):
        return
    sc.pp.highly_variable_genes(
        adata, flavor="seurat_v3", n_top_genes=n_top_genes,
        layer="counts",
        batch_key=batch_key if (batch_key and batch_key in adata.obs) else None,
    )


def train_scvi(
    input_h5ad: str = "calcinosis_qc.h5ad",
    out_basename: str = "calcinosis_scvi",
    n_top_hvgs: int = N_TOP_HVGS,
    batch_key: str | None = "sample",
    n_latent: int = N_LATENT,
    max_epochs: int = MAX_EPOCHS,
) -> sc.AnnData:

    # ── seeds ──────────────────────────────────────────────────────────────
    np.random.seed(SEEDS["numpy"])
    scvi.settings.seed = SEEDS["scvi"]

    # ── load ───────────────────────────────────────────────────────────────
    adata = sc.read_h5ad(PROC / input_h5ad)
    print(f"[scVI] Loaded {adata.n_obs} cells × {adata.n_vars} genes")

    # Ensure counts layer
    if "counts" not in adata.layers:
        raise ValueError("No 'counts' layer found. Re-run preprocess.py.")

    # HVGs
    _ensure_hvgs(adata, n_top_hvgs, batch_key)
    adata_hvg = adata[:, adata.var["highly_variable"]].copy()
    print(f"[scVI] Using {adata_hvg.n_vars} HVGs for training")

    # ── scVI setup & training ──────────────────────────────────────────────
    scvi.model.SCVI.setup_anndata(
        adata_hvg,
        layer="counts",
        batch_key=batch_key if (batch_key and batch_key in adata_hvg.obs) else None,
    )
    model = scvi.model.SCVI(
        adata_hvg,
        n_latent=n_latent,
        n_layers=N_LAYERS,
        n_hidden=N_HIDDEN,
    )
    print("[scVI] Training...")
    model.train(
        max_epochs=max_epochs,
        batch_size=BATCH_SIZE,
        early_stopping=EARLY_STOP,
        plan_kwargs={"lr": 1e-3},
    )

    # ── embed & cluster ───────────────────────────────────────────────────
    adata_hvg.obsm["X_scVI"] = model.get_latent_representation()
    sc.pp.neighbors(adata_hvg, use_rep="X_scVI", n_neighbors=15)
    sc.tl.umap(adata_hvg, random_state=SEEDS["numpy"])
    sc.tl.leiden(adata_hvg, resolution=LEIDEN_RES, random_state=SEEDS["numpy"])
    print(f"[scVI] {adata_hvg.obs['leiden'].nunique()} Leiden clusters")

    # ── transfer embeddings to full adata ─────────────────────────────────
    adata.obsm["X_scVI"] = adata_hvg.obsm["X_scVI"]
    adata.obsm["X_umap"] = adata_hvg.obsm["X_umap"]
    adata.obs["leiden"]  = adata_hvg.obs["leiden"]
    adata.obsp           = adata_hvg.obsp
    adata.uns            = adata_hvg.uns

    # ── calcinosis gene signatures ─────────────────────────────────────────
    symbols = _symbol_series(adata)
    var_symbols_upper = symbols.str.upper().values

    for sig_name, sig_genes in CALCINOSIS_SIGNATURES.items():
        present = [g for g in sig_genes
                   if g.upper() in var_symbols_upper]
        if present:
            sc.tl.score_genes(adata, present,
                              score_name=f"score_{sig_name.lower()}")
            print(f"  Scored {sig_name}: {len(present)}/{len(sig_genes)} genes found")
        else:
            print(f"  [warn] No genes found for {sig_name} signature")

    # ── save ──────────────────────────────────────────────────────────────
    out = PROC / f"{out_basename}.h5ad"
    adata.write_h5ad(out)
    print(f"[saved] {out}")

    # Save model
    model_dir = PROC / f"{out_basename}_model"
    model.save(str(model_dir), overwrite=True)
    print(f"[saved] scVI model -> {model_dir}")

    return adata


def annotate(
    input_h5ad: str = "calcinosis_scvi.h5ad",
    out_basename: str = "calcinosis_scvi_annot",
) -> sc.AnnData:
    """Score marker genes and assign coarse cell-type labels."""

    adata = sc.read_h5ad(PROC / input_h5ad)
    symbols = _symbol_series(adata)
    var_symbols_upper = symbols.str.upper().values

    scores = {}
    for ct, genes in MARKERS.items():
        present = [g for g in genes if g.upper() in var_symbols_upper]
        if present:
            sc.tl.score_genes(adata, present, score_name=f"_ct_{ct}")
            scores[ct] = adata.obs[f"_ct_{ct}"].values
        else:
            scores[ct] = np.zeros(adata.n_obs)

    score_df = pd.DataFrame(scores, index=adata.obs_names)
    adata.obs["cell_type"] = score_df.idxmax(axis=1)

    # Per-cluster majority vote
    cluster_ct = (
        adata.obs.groupby("leiden")["cell_type"]
        .agg(lambda x: x.value_counts().index[0])
    )
    adata.obs["cell_type_cluster"] = adata.obs["leiden"].map(cluster_ct)

    # UMAP figures
    for colour in ["leiden", "condition", "cell_type_cluster",
                   "score_osteogenic", "score_calcification",
                   "score_hypoxia", "score_fibrosis"]:
        if colour in adata.obs.columns or colour in adata.obsm:
            try:
                sc.pl.umap(adata, color=colour, show=False,
                           save=f"_{colour}.png")
            except Exception:
                pass

    out = PROC / f"{out_basename}.h5ad"
    adata.write_h5ad(out)
    print(f"[saved] {out}")
    return adata


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input",  default="calcinosis_qc.h5ad")
    ap.add_argument("--out",    default="calcinosis_scvi")
    ap.add_argument("--epochs", type=int, default=MAX_EPOCHS)
    args = ap.parse_args()

    print("=" * 60)
    print("Calcinosis / CREST Paper 01 — scVI Modelling")
    print("=" * 60)

    adata = train_scvi(args.input, args.out, max_epochs=args.epochs)
    annotate(f"{args.out}.h5ad", f"{args.out}_annot")


if __name__ == "__main__":
    main()
