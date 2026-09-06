"""Raw-count QC, sample-level doublet detection and descriptive embeddings."""

from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc

from .core import audit, mkdir, scores, validate_adata


def run(cfg):
    path = Path(cfg["input_dir"]) / "singlecell.h5ad"
    out = mkdir(Path(cfg["output_dir"]) / "singlecell")
    a = ad.read_h5ad(path)
    validate_adata(a, "singlecell", cfg["mode"])
    if not cfg["analysis"]["allow_unreviewed_annotations"]:
        rejected = {"unknown", "unreviewed", "synthetic_truth", "automatic", ""}
        if a.obs.annotation_source.astype(str).str.lower().isin(rejected).any():
            raise ValueError("Real analysis requires reviewed cell labels and annotation_source")
    a.X = a.layers["counts"].copy()
    a.var["mt"] = a.var.gene_symbol.str.upper().str.startswith("MT-").to_numpy()
    sc.pp.calculate_qc_metrics(a, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)
    qc = cfg["qc"]
    a.obs["pass_qc"] = ((a.obs.n_genes_by_counts >= qc["min_genes"])
                        & (a.obs.n_genes_by_counts <= qc["max_genes"])
                        & (a.obs.pct_counts_mt <= qc["max_mt_pct"]))
    a.obs["doublet_score"] = np.nan
    a.obs["predicted_doublet"] = False
    if qc["doublets"]:
        # Fit separately to each capture; mixing donors changes the simulated null.
        for _, obs in a.obs.groupby("sample_id", observed=True):
            b = a[obs.index[a.obs.loc[obs.index, "pass_qc"]]].copy()
            if b.n_obs < qc["min_cells_per_sample"]:
                raise ValueError("Too few cells for a capture: inspect exclusions before rerunning")
            sc.pp.scrublet(b, expected_doublet_rate=qc.get("expected_doublet_rate", .06),
                          random_state=cfg["seed"], verbose=False)
            a.obs.loc[b.obs_names, "doublet_score"] = b.obs.doublet_score
            a.obs.loc[b.obs_names, "predicted_doublet"] = b.obs.predicted_doublet
    a.obs["retained"] = a.obs.pass_qc & ~a.obs.predicted_doublet
    a.obs.to_csv(out / "cell_qc.csv")
    original_samples = set(a.obs.sample_id)
    a = a[a.obs.retained].copy()
    if set(a.obs.sample_id) != original_samples:
        raise ValueError("A complete capture was lost to QC; review the cohort and exclusions")
    counts = a.obs.groupby("sample_id", observed=True).size()
    if counts.empty or counts.min() < qc["min_cells_per_sample"]:
        raise ValueError("Capture below min_cells_per_sample after QC; inspect cell_qc.csv")
    sc.pp.filter_genes(a, min_cells=qc["min_cells_per_gene"])
    a.layers["counts"] = a.X.copy()
    sc.pp.normalize_total(a, target_sum=1e4)
    sc.pp.log1p(a)
    sc.pp.highly_variable_genes(a, flavor="seurat", n_top_genes=min(
        cfg["analysis"]["hvg"], a.n_vars - 1), batch_key="sample_id")
    n_hvg = int(a.var.highly_variable.sum())
    n_pcs = min(cfg["analysis"]["pcs"], n_hvg - 1, a.n_obs - 1)
    if n_pcs < 2:
        raise ValueError("Insufficient variable genes/cells for PCA")
    # Counts are untouched. Embedding alone is never used for patient-level DE.
    sc.pp.pca(a, n_comps=n_pcs, mask_var="highly_variable", random_state=cfg["seed"])
    sc.pp.neighbors(a, n_neighbors=min(cfg["analysis"]["neighbors"], a.n_obs - 1),
                    n_pcs=n_pcs, random_state=cfg["seed"])
    sc.tl.umap(a, random_state=cfg["seed"])
    sc.tl.leiden(a, resolution=cfg["analysis"]["resolution"], random_state=cfg["seed"],
                 flavor="leidenalg", n_iterations=2, directed=False)
    program, coverage = scores(a, cfg["analysis"]["min_signature_coverage"])
    for name in program:
        a.obs[f"score_{name}"] = program[name]
    coverage.to_csv(out / "program_coverage.csv", index=False)
    composition = pd.crosstab(a.obs.patient_id, a.obs.cell_type)
    composition.div(composition.sum(axis=1), axis=0).to_csv(out / "patient_composition.csv")
    a.obs.groupby(["patient_id", "condition", "cell_type"], observed=True)[
        [f"score_{p}" for p in program]].mean().to_csv(out / "patient_programs.csv")
    a.write_h5ad(out / "processed.h5ad", compression="gzip")
    a.var.to_csv(out / "genes.csv")
    audit(cfg, "singlecell", [path], {"cells_retained": a.n_obs,
          "patients": a.obs.patient_id.nunique(), "genes": a.n_vars,
          "embedding": "uncorrected PCA; inspect batch and biology before integration"})
