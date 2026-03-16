import scanpy as sc
import scvi
adata = sc.read_h5ad("calcinosis-crest/data/processed/calcinosis_qc.h5ad")
sc.pp.highly_variable_genes(adata, n_top_genes=4000, flavor="seurat_v3", layer="counts", batch_key="sample")
adata_hvg = adata[:, adata.var["highly_variable"]].copy()
scvi.model.SCVI.setup_anndata(adata_hvg, layer="counts", batch_key="sample")
model = scvi.model.SCVI(adata_hvg, n_latent=30, n_layers=2, n_hidden=128)
model.train(max_epochs=2, batch_size=256, early_stopping=False, plan_kwargs={"lr": 1e-3})
