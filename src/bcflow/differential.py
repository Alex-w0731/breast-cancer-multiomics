"""Patient-pseudobulk and bulk negative-binomial differential expression."""

from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from formulaic import model_matrix
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
from scipy import sparse

from .core import audit, bh, mkdir, read_count_csv, validate_counts, write_json


def aggregate(a, min_cells):
    """Sum raw UMIs per patient/cell type, including technical replicates once.

    A patient spanning multiple batches is not silently assigned a batch. Such
    data need an explicitly modeled repeated-measure design outside this engine.
    """
    meta, rows, dropped = [], [], []
    for (patient, ct), frame in a.obs.groupby(["patient_id", "cell_type"], observed=True):
        if len(frame) < min_cells:
            dropped.append({"patient_id": patient, "cell_type": ct, "n_cells": len(frame)})
            continue
        for col in ["condition", "batch"]:
            if frame[col].nunique() != 1:
                raise ValueError(f"Patient {patient} has multiple {col} values; model explicitly")
        idx = a.obs_names.get_indexer(frame.index)
        row = sparse.csr_matrix(a.layers["counts"])[idx].sum(axis=0)
        rows.append(np.asarray(row).ravel().astype(np.int64))
        record = {"sample_id": f"{patient}::{ct}", "patient_id": str(patient),
                  "cell_type": str(ct), "n_cells": len(frame)}
        # Preserve donor-constant covariates, permitting an audited adjusted design.
        for col in frame:
            if col not in record and frame[col].nunique(dropna=False) == 1:
                record[col] = frame[col].iloc[0]
        record["sample_id"] = f"{patient}::{ct}"
        meta.append(record)
    if not rows:
        raise ValueError("No eligible pseudobulk samples")
    md = pd.DataFrame(meta).set_index("sample_id")
    return pd.DataFrame(rows, index=md.index, columns=a.var_names), md, pd.DataFrame(dropped)


def check_design(counts, meta, cfg):
    key, test, control = cfg["contrast"]
    for col in ["patient_id", key]:
        if col not in meta or meta[col].isna().any() or meta[col].astype(str).str.strip().eq("").any():
            raise ValueError(f"Missing or blank metadata: {col}")
    if not counts.index.is_unique or not counts.columns.is_unique or not meta.index.is_unique:
        raise ValueError("Duplicate sample/gene IDs")
    if not counts.index.equals(meta.index):
        raise ValueError("Counts and metadata must have identical, ordered sample IDs")
    if meta.patient_id.duplicated().any():
        raise ValueError("One independent observation per patient is required for this model")
    validate_counts(counts.to_numpy())
    if set(meta[key].astype(str)) != {test, control}:
        raise ValueError("Exactly the two prespecified contrast levels must be present")
    n = meta.groupby(key, observed=True).patient_id.nunique()
    if n.min() < cfg["analysis"]["min_donors_per_group"]:
        raise ValueError("Insufficient independent patients per contrast group")
    design = model_matrix(cfg["design"], meta, na_action="raise")
    design = np.asarray(design, dtype=float)
    if not np.isfinite(design).all() or np.linalg.matrix_rank(design) < design.shape[1]:
        raise ValueError("Rank-deficient design: condition/batch confounding or redundant covariates")
    if design.shape[0] - design.shape[1] < 2:
        raise ValueError("Insufficient residual degrees of freedom")
    # Ensure the contrast factor is represented in the fitted model.
    if key not in cfg["design"].replace("~", " ").replace("+", " ").split():
        raise ValueError("Contrast factor missing from design")


def fit(counts, meta, cfg):
    check_design(counts, meta, cfg)
    spec = cfg["analysis"]
    keep = (counts >= spec["min_count"]).sum(axis=0) >= spec["min_samples_gene"]
    counts = counts.loc[:, keep].astype(np.int64)
    if counts.shape[1] < 20:
        raise ValueError("Fewer than 20 genes survive count filtering")
    dds = DeseqDataSet(counts=counts, metadata=meta, design=cfg["design"],
                       refit_cooks=True, n_cpus=spec["threads"], quiet=True)
    dds.deseq2()
    ds = DeseqStats(dds, contrast=cfg["contrast"], n_cpus=spec["threads"], quiet=True,
                    cooks_filter=True, independent_filter=True)
    ds.summary()
    result = ds.results_df.copy()
    result.index.name = "gene_id"
    # Raw Wald coefficients, not mislabeled as shrinkage estimates.
    result["estimate_type"] = "unshrunk_Wald_log2FC"
    normalized = pd.DataFrame(dds.layers["normed_counts"], index=counts.index,
                               columns=counts.columns)
    return result, normalized


def pseudobulk(cfg):
    path = Path(cfg["output_dir"]) / "singlecell/processed.h5ad"
    out = mkdir(Path(cfg["output_dir"]) / "pseudobulk")
    a = ad.read_h5ad(path)
    counts, md, dropped = aggregate(a, cfg["analysis"]["min_cells_pseudobulk"])
    counts.to_csv(out / "counts.csv")
    md.to_csv(out / "metadata.csv")
    dropped.to_csv(out / "excluded_low_cell_count.csv", index=False)
    results, statuses = [], []
    for ct in sorted(md.cell_type.unique()):
        sub = md[md.cell_type == ct].copy()
        n = sub.groupby("condition", observed=True).patient_id.nunique()
        if len(n) != 2 or n.min() < cfg["analysis"]["min_donors_per_group"]:
            statuses.append({"cell_type": ct, "status": "not_testable", "donors": n.to_dict()})
            continue
        res, _ = fit(counts.loc[sub.index], sub, cfg)
        res["cell_type"] = ct
        res["gene_symbol"] = a.var.gene_symbol.reindex(res.index)
        results.append(res.reset_index())
        statuses.append({"cell_type": ct, "status": "tested", "donors": n.to_dict()})
    write_json(out / "test_status.json", statuses)
    if not results:
        raise ValueError("No cell type meets independent-donor threshold; inspect test_status.json")
    result = pd.concat(results, ignore_index=True)
    # Primary multiple-testing family: all gene x cell-type contrasts together.
    result["q_global"] = bh(result.pvalue)
    result.to_csv(out / "differential.csv", index=False)
    audit(cfg, "pseudobulk", [path], {"test_status": statuses, "fdr_family": "all gene x celltype"})


def bulk(cfg):
    inp = Path(cfg["input_dir"])
    out = mkdir(Path(cfg["output_dir"]) / "bulk")
    counts = read_count_csv(inp / "bulk_counts.csv")
    md = pd.read_csv(inp / "bulk_metadata.csv", index_col=0)
    genes = pd.read_csv(inp / "bulk_genes.csv", index_col=0)
    if cfg["mode"] == "real" and "data_status" in md:
        if md.data_status.eq("SYNTHETIC").any():
            raise ValueError("Synthetic bulk records cannot enter a real run")
    if not genes.index.is_unique or not set(counts.columns).issubset(genes.index):
        raise ValueError("Missing or ambiguous bulk gene annotation")
    if set(counts.index) != set(md.index):
        raise ValueError("Bulk counts and metadata sample sets differ")
    md = md.loc[counts.index]
    result, normalized = fit(counts, md, cfg)
    result["gene_symbol"] = genes.gene_symbol.reindex(result.index)
    result.to_csv(out / "differential.csv")
    normalized.to_csv(out / "normalized_counts.csv")
    md.to_csv(out / "metadata.csv")
    audit(cfg, "bulk", [inp / x for x in ["bulk_counts.csv", "bulk_metadata.csv", "bulk_genes.csv"]],
          {"patients": md.patient_id.nunique(), "design": cfg["design"],
           "fdr_family": "bulk genes", "normalization": "DESeq2 size factors"})
