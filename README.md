# Calcinosis in Systemic Sclerosis — Single-Cell Atlas and Drug Repurposing
### CREST Implication Testing Series — Paper 01 (C = Calcinosis)

**DOI:** *(assigned after Zenodo upload)*

This repository contains a reproducible single-cell RNA-seq pipeline to
identify the cellular and molecular mechanisms of calcinosis in systemic
sclerosis (SSc), and to prioritise drug repurposing candidates via LINCS
L1000 transcriptomic reversal scoring.

Calcinosis — the deposition of calcium hydroxyapatite in skin and soft
tissue — affects 20–40% of SSc patients and has no approved treatment.
This paper applies the same scVI + Wilcoxon DE + Enrichr/LINCS pipeline
used in our fibrosis paper (Ritschel & Claude, 2026) to the calcinosis
problem using the higher-resolution GSE195452 fibroblast atlas.

---

## Dataset

**GEO accession:** [GSE195452](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE195452)

Gur C et al. (2022). *LGR5 expressing skin fibroblasts define a major cellular
hub perturbed in scleroderma.* **Cell**, 185(8), 1373–1388.
https://doi.org/10.1016/j.cell.2022.03.011

GSE195452 is a single-cell multiome (RNA + ATAC) atlas of SSc and healthy
skin with 49 fibroblast subsets. We use the RNA modality. This dataset was
chosen because:
- It resolves LGR5+ fibroblast progenitors — the cell type most likely to
  undergo osteogenic differentiation leading to calcinosis
- It includes both healthy and SSc samples enabling SSc vs Healthy DE
- It has been widely cited in calcinosis-focused SSc papers
- It is publicly available without institutional access requirements

---

## Pipeline

```
python3 src/data_download.py    # download GSE195452 from GEO
python3 src/preprocess.py       # QC, normalise, HVG selection
python3 src/modeling.py         # scVI embedding, UMAP, Leiden, calcinosis scoring
python3 src/de_analysis.py      # Wilcoxon DE + Enrichr + LINCS L1000
```

### Key analysis steps

| Step | Script | Output |
|------|--------|--------|
| Download | `data_download.py` | `data/raw/GSE195452/*.h5` |
| QC + preprocessing | `preprocess.py` | `calcinosis_qc.h5ad` |
| scVI + clustering | `modeling.py` | `calcinosis_scvi_annot.h5ad` |
| Calcinosis signatures | `modeling.py` | scores: osteogenic, calcification, hypoxia |
| Wilcoxon DE | `de_analysis.py` | `de_leiden_wilcoxon.csv` |
| SSc vs Healthy DE | `de_analysis.py` | `de_ssc_vs_healthy_per_cluster.csv` |
| Enrichr (GO/Reactome/KEGG) | `de_analysis.py` | `results/drug_repurposing/leiden_*/` |
| LINCS L1000 reversal | `de_analysis.py` | `lincs_reversal_top15_by_cluster.csv` |

---

## Calcinosis Biology Focus

The analysis specifically scores and tracks genes relevant to calcinosis
pathophysiology:

**Osteogenic differentiation:** RUNX2, SP7, BGLAP, ALPL, COL10A1, IBSP

**Calcification pathway:** ENPP1, ANKH, MGP, FETUB, BMP2/4, SMAD1/5

**Hypoxia response:** HIF1A, VEGFA, LDHA, SLC2A1, BNIP3
(hypoxia is a primary driver of calcinosis; Valenzuela et al. 2016)

**Fibrosis:** TGFB1/2, CTGF, FN1, POSTN, THBS1

**Inflammation:** IL6, IL1B, TNF, CXCL8, CXCL12, CCL2

---

## Environment

```bash
# Same environment as paper 01 (scleroderma-scvi)
mamba env create -f environment.yml
conda activate ssc-scvi
```

Key versions:
- Python: 3.10
- scanpy: 1.11.4
- scvi-tools: 1.3.3
- gseapy: 1.1.10
- anndata: 0.11.4

---

## LINCS Libraries Used

| Library | Purpose |
|---------|---------|
| GO_Biological_Process_2023 | Pathway enrichment |
| Reactome_2022 | Pathway enrichment |
| KEGG_2021_Human | Pathway enrichment |
| LINCS_L1000_Chem_Pert_up | Drug repurposing (UP perturbation) |
| LINCS_L1000_Chem_Pert_down | Drug repurposing (DOWN perturbation) |

Reversal score = sign × (−log₁₀ adjusted p-value), where sign is +1 when
the drug perturbation direction opposes the disease gene expression direction.

---

## Key References

- Gur C et al. (2022). LGR5 expressing skin fibroblasts define a major
  cellular hub perturbed in scleroderma. *Cell*, 185, 1373–1388.
- Valenzuela A et al. (2016). Calcinosis is associated with digital ulcers
  and osteoporosis in SSc. *Semin Arthritis Rheum*, 46, 344–349.
- Herrick AL & Gallas A. (2016). Systemic sclerosis-related calcinosis.
  *J Scleroderma Relat Disord*, 1, 194–203.
- Ritschel G & Claude (Anthropic). (2026). Single-cell atlas of systemic
  sclerosis skin reveals therapeutic targets. Zenodo.
  https://doi.org/10.5281/zenodo.17156327

---

## Citation

```bibtex
@misc{ritschel_claude_2026_calcinosis_scrna,
  author    = {Ritschel, Glen and Claude (Anthropic)},
  title     = {Single-cell transcriptomic analysis of calcinosis in
               systemic sclerosis identifies osteogenic fibroblast states
               and drug repurposing candidates},
  year      = {2026},
  publisher = {Zenodo},
  doi       = {10.5281/zenodo.XXXXXXX}
}
```

---

## CREST Series

| Paper | Letter | Manifestation | Status |
|-------|--------|--------------|--------|
| 01 | **C** | Calcinosis | **This paper** |
| 02 | **R** | Raynaud phenomenon | Planned |
| 03 | **E** | Esophageal dysmotility | Planned |
| 04 | **S** | Sclerodactyly | Planned |
| 05 | **T** | Telangiectasia | Planned |
