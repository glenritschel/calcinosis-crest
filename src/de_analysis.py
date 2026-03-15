#!/usr/bin/env python3
"""
de_analysis.py  — Calcinosis / CREST Paper 01
----------------------------------------------
Wilcoxon differential expression + Enrichr (GSEA) + LINCS L1000
drug-repurposing analysis for GSE195452 calcinosis paper.

Two DE comparisons run in parallel:
    A. Per-Leiden-cluster (same as paper 01)
    B. SSc vs Healthy within each cluster (calcinosis-specific angle)

Then for each cluster:
    - Top 200 UP / DOWN genes → Enrichr (GO_BP, Reactome, KEGG, LINCS L1000)
    - Signed reversal score = -log10(adj_p) * sign  (LINCS drug repurposing)

Usage:
    python3 src/de_analysis.py

Outputs:
    results/tables/de_leiden_wilcoxon.csv
    results/tables/de_ssc_vs_healthy_per_cluster.csv
    results/drug_repurposing/leiden_*/up/enrichr_*.csv
    results/drug_repurposing/leiden_*/down/enrichr_*.csv
    results/drug_repurposing/lincs_reversal_top15_by_cluster.csv

Environment: ssc-scvi  (gseapy 1.1.10, scanpy 1.11.4)
"""

from pathlib import Path
import math
import re
import time
import numpy as np
import pandas as pd
import scanpy as sc
import gseapy as gp

PROC  = Path("data/processed")
RES   = Path("results/tables");           RES.mkdir(parents=True, exist_ok=True)
DR    = Path("results/drug_repurposing"); DR.mkdir(parents=True, exist_ok=True)

# Enrichr libraries — identical to paper 01
ENRICHR_LIBS = [
    "GO_Biological_Process_2023",
    "Reactome_2022",
    "KEGG_2021_Human",
    "LINCS_L1000_Chem_Pert_up",
    "LINCS_L1000_Chem_Pert_down",
]

N_TOP_GENES = 200   # genes per direction submitted to Enrichr
MIN_CELLS   = 20    # minimum cells in a cluster to run DE


# ── helpers ──────────────────────────────────────────────────────────────────

def _symbol_series(adata: sc.AnnData) -> pd.Series:
    for c in ["gene_symbol", "gene_symbols", "symbol", "SYMBOL",
              "gene_name", "gene_names"]:
        if c in adata.var.columns:
            return adata.var[c].astype(str)
    return adata.var_names.to_series().astype(str)


def _to_symbols(gene_list: list, sym_map: dict) -> list:
    """Convert var_names (possibly Ensembl) to symbols where possible."""
    return [sym_map.get(g, g) for g in gene_list]


def _enrichr_safe(gene_list: list, lib: str, description: str,
                  outdir: Path, retries: int = 3) -> pd.DataFrame | None:
    """Run gseapy.enrichr with retries; return results df or None."""
    if not gene_list:
        return None
    outdir.mkdir(parents=True, exist_ok=True)
    for attempt in range(retries):
        try:
            enr = gp.enrichr(
                gene_list=gene_list,
                gene_sets=lib,
                organism="human",
                
                outdir=str(outdir),
                verbose=False,
            )
            return enr.results
        except Exception as e:
            print(f"    [warn] Enrichr attempt {attempt+1}/{retries} failed "
                  f"for {lib}: {e}")
            time.sleep(3 * (attempt + 1))
    return None


def _signed_reversal_score(row, src: str) -> float:
    """
    + if cluster UP genes overlap LINCS drug DOWN perturbation  (reversal)
    + if cluster DOWN genes overlap LINCS drug UP perturbation  (reversal)
    - opposite pairings
    Score magnitude = -log10(adj_p)
    """
    lib  = str(row.get("Gene_set", row.get("library", ""))).lower()
    adjp = row.get("Adjusted P-value", row.get("adj_p", 1.0))
    try:
        adjp = float(adjp)
    except Exception:
        adjp = 1.0
    if adjp <= 0 or np.isnan(adjp):
        adjp = 1.0
    score = -math.log10(adjp)
    sign  = 0
    if src == "up":
        sign = +1 if "chem_pert_down" in lib else (-1 if "chem_pert_up" in lib else 0)
    elif src == "down":
        sign = +1 if "chem_pert_up" in lib else (-1 if "chem_pert_down" in lib else 0)
    return sign * score


def _clean_compound(term: str) -> str:
    t = str(term).replace("|", " ").replace(";", " ")
    t = re.split(r"[_\[]", t)[0]
    t = re.sub(r"\s*\(.*?\)\s*$", "", t)
    return t.strip()


# ── DE per Leiden cluster ─────────────────────────────────────────────────────

def de_per_cluster(adata: sc.AnnData, sym_map: dict) -> pd.DataFrame:
    """Wilcoxon DE — each cluster vs rest."""
    sc.tl.rank_genes_groups(adata, groupby="leiden", method="wilcoxon",
                             use_raw=False)
    df = sc.get.rank_genes_groups_df(adata, group=None)
    # Attach symbols
    df["gene_symbol"] = df["names"].map(sym_map).fillna(df["names"])
    out = RES / "de_leiden_wilcoxon.csv"
    df.to_csv(out, index=False)
    print(f"[saved] {out}  ({len(df)} rows)")
    return df


# ── DE SSc vs Healthy per cluster ─────────────────────────────────────────────

def de_ssc_vs_healthy(adata: sc.AnnData, sym_map: dict) -> pd.DataFrame:
    """Wilcoxon DE between SSc and Healthy within each Leiden cluster."""
    if "condition" not in adata.obs.columns:
        print("[warn] No 'condition' column — skipping SSc vs Healthy DE")
        return pd.DataFrame()

    results = []
    for clus in sorted(adata.obs["leiden"].unique()):
        mask = adata.obs["leiden"] == clus
        sub  = adata[mask]
        n_ssc     = (sub.obs["condition"] == "SSc").sum()
        n_healthy = (sub.obs["condition"] == "Healthy").sum()
        if n_ssc < MIN_CELLS or n_healthy < MIN_CELLS:
            print(f"  [skip cluster {clus}] SSc={n_ssc}, Healthy={n_healthy}")
            continue
        sc.tl.rank_genes_groups(sub, groupby="condition", groups=["SSc"],
                                 reference="Healthy", method="wilcoxon",
                                 use_raw=False)
        df = sc.get.rank_genes_groups_df(sub, group="SSc")
        df["leiden"]       = clus
        df["gene_symbol"]  = df["names"].map(sym_map).fillna(df["names"])
        results.append(df)

    if not results:
        return pd.DataFrame()

    out_df = pd.concat(results, ignore_index=True)
    out    = RES / "de_ssc_vs_healthy_per_cluster.csv"
    out_df.to_csv(out, index=False)
    print(f"[saved] {out}  ({len(out_df)} rows)")
    return out_df


# ── Enrichr + LINCS per cluster ───────────────────────────────────────────────

def run_enrichr_lincs(de_df: pd.DataFrame,
                      adata: sc.AnnData,
                      sym_map: dict) -> pd.DataFrame:
    """
    For each Leiden cluster, submit top UP/DOWN gene symbols to Enrichr
    across all ENRICHR_LIBS. Compute LINCS reversal scores.
    """
    ct_map = {}
    if "cell_type_cluster" in adata.obs.columns:
        ct_map = adata.obs.groupby("leiden")["cell_type_cluster"].first().to_dict()

    records = []
    clusters = sorted(de_df["group"].unique() if "group" in de_df.columns
                      else de_df["leiden"].unique()
                      if "leiden" in de_df.columns else [])

    for clus in clusters:
        if "group" in de_df.columns:
            sub = de_df[de_df["group"] == str(clus)]
        else:
            sub = de_df[de_df["leiden"] == str(clus)]

        # Sort by score (logfoldchange or scores column)
        score_col = next((c for c in ["scores", "logfoldchanges", "score"]
                          if c in sub.columns), None)
        if score_col:
            sub = sub.sort_values(score_col, ascending=False)

        sym_col = "gene_symbol" if "gene_symbol" in sub.columns else "names"
        up_genes   = sub.head(N_TOP_GENES)[sym_col].tolist()
        down_genes = sub.tail(N_TOP_GENES)[sym_col].tolist()

        cell_type = ct_map.get(str(clus), f"cluster_{clus}")
        clus_dir  = DR / f"leiden_{clus}"

        for direction, gene_list in [("up", up_genes), ("down", down_genes)]:
            gene_list = [g for g in gene_list if isinstance(g, str) and g]
            if not gene_list:
                continue
            out_sub = clus_dir / direction
            for lib in ENRICHR_LIBS:
                desc = f"leiden{clus}_{direction}"
                print(f"  Enrichr cluster={clus} dir={direction} lib={lib[:30]}...")
                res  = _enrichr_safe(gene_list, lib, desc, out_sub)
                if res is None or res.empty:
                    continue
                # Save per-library CSV
                csv = out_sub / f"enrichr_leiden{clus}_{direction}_{lib[:20]}.csv"
                res.to_csv(csv, index=False)

                # Collect LINCS reversal records
                if "LINCS" in lib:
                    for _, row in res.iterrows():
                        score = _signed_reversal_score(row, direction)
                        if score == 0:
                            continue
                        records.append({
                            "cluster":    str(clus),
                            "cell_type":  cell_type,
                            "direction":  direction,
                            "library":    lib,
                            "term":       row.get("Term", ""),
                            "compound":   _clean_compound(row.get("Term", "")),
                            "adj_p":      row.get("Adjusted P-value", 1.0),
                            "rev_score":  score,
                        })

    if not records:
        print("[warn] No LINCS records collected.")
        return pd.DataFrame()

    df_lincs = pd.DataFrame(records)

    # Aggregate: top 15 per cluster by max reversal score
    agg = (
        df_lincs.groupby(["cluster", "cell_type", "compound"])
        .agg(
            rev_score_sum=("rev_score", "sum"),
            rev_score_max=("rev_score", "max"),
            best_p       =("adj_p",     "min"),
            hits         =("rev_score", "count"),
        )
        .reset_index()
    )
    top15 = (
        agg.sort_values("rev_score_max", ascending=False)
        .groupby("cluster")
        .head(15)
        .reset_index(drop=True)
    )
    out = DR / "lincs_reversal_top15_by_cluster.csv"
    top15.to_csv(out, index=False)
    print(f"[saved] {out}  ({len(top15)} rows)")

    # Cross-cluster aggregate
    cross = (
        agg.groupby("compound")
        .agg(
            n_clusters   =("cluster",       "nunique"),
            total_score  =("rev_score_sum",  "sum"),
            max_rev      =("rev_score_max",  "max"),
            best_p       =("best_p",         "min"),
        )
        .reset_index()
        .sort_values("total_score", ascending=False)
    )
    out2 = DR / "lincs_reversal_crosscluster.csv"
    cross.to_csv(out2, index=False)
    print(f"[saved] {out2}")

    return top15


# ── main ──────────────────────────────────────────────────────────────────────

def main():
    print("=" * 60)
    print("Calcinosis / CREST Paper 01 — DE + GSEA + LINCS")
    print("=" * 60)

    # Load latest annotated h5ad
    cands = sorted(
        PROC.glob("calcinosis_scvi*annot*.h5ad"),
        key=lambda p: p.stat().st_mtime, reverse=True
    )
    if not cands:
        cands = sorted(
            PROC.glob("calcinosis_scvi*.h5ad"),
            key=lambda p: p.stat().st_mtime, reverse=True
        )
    if not cands:
        raise FileNotFoundError(
            "No processed h5ad found. Run preprocess.py then modeling.py first."
        )
    adata_path = cands[0]
    print(f"[info] Using AnnData: {adata_path}")
    adata = sc.read_h5ad(adata_path)
    print(f"[info] {adata.n_obs} cells, {adata.n_vars} genes, "
          f"{adata.obs['leiden'].nunique()} clusters")

    # Symbol map  (var_names -> gene_symbol)
    sym_series = _symbol_series(adata)
    sym_map    = dict(zip(adata.var_names, sym_series))

    # 1. DE per cluster
    print("\n[1] Wilcoxon DE per Leiden cluster...")
    de_cluster = de_per_cluster(adata, sym_map)

    # 2. DE SSc vs Healthy (calcinosis angle)
    print("\n[2] Wilcoxon DE SSc vs Healthy per cluster...")
    de_ssc = de_ssc_vs_healthy(adata, sym_map)

    # 3. Enrichr + LINCS on per-cluster DE
    print("\n[3] Enrichr + LINCS L1000 drug repurposing...")
    # Use SSc-vs-Healthy DE if available; fall back to cluster DE
    de_for_lincs = de_ssc if (not de_ssc.empty and "leiden" in de_ssc.columns) \
                   else de_cluster
    top15 = run_enrichr_lincs(de_for_lincs, adata, sym_map)

    print("\n[done] de_analysis.py complete.")
    print(f"  DE rows:    {len(de_cluster)}")
    print(f"  SSc DE rows:{len(de_ssc)}")
    print(f"  LINCS top15:{len(top15)}")


if __name__ == "__main__":
    main()
