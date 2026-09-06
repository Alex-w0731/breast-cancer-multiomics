"""Section-aware spatial QC, reference mapping and patient-level summaries."""

from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.optimize import nnls
from sklearn.neighbors import NearestNeighbors

from .core import audit, bh, dense, mkdir, scores, validate_adata, write_json


def graph(coords, k=6):
    """Symmetric local graph for ONE section; avoid long bridges across empty tissue."""
    coords = np.asarray(coords, dtype=float)
    if len(coords) < 4 or len(np.unique(coords, axis=0)) != len(coords):
        raise ValueError("Section needs >=4 distinct spatial positions")
    distances, idx = NearestNeighbors(n_neighbors=min(k + 1, len(coords))).fit(
        coords).kneighbors(coords)
    radius = 1.8 * np.median(distances[:, 1])
    rows = np.repeat(np.arange(len(coords)), idx.shape[1] - 1)
    cols = idx[:, 1:].ravel()
    keep = distances[:, 1:].ravel() <= radius
    w = sparse.csr_matrix((np.ones(keep.sum()), (rows[keep], cols[keep])),
                           shape=(len(coords), len(coords)))
    w = w.maximum(w.T)
    w.setdiag(0)
    w.eliminate_zeros()
    if w.nnz == 0:
        raise ValueError("No spatial neighbors under the gap-aware graph")
    return w


def moran(values, w, permutations, rng):
    """Two-sided permutation test against spatial randomness, within one section."""
    z = np.asarray(values, dtype=float) - np.mean(values)
    denominator = z @ z
    if denominator <= 0:
        return np.nan, np.nan
    factor = len(z) / w.sum()
    observed = float(factor * (z @ (w @ z)) / denominator)
    expected = -1 / (len(z) - 1)
    null = []
    for _ in range(permutations):
        shuffled = rng.permutation(z)
        null.append(float(factor * (shuffled @ (w @ shuffled)) / denominator))
    p = (1 + np.sum(np.abs(np.asarray(null) - expected) >= abs(observed - expected))) / (
        permutations + 1)
    return observed, float(p)


def reference_signature(ref, genes):
    """Average CPM by patient within type, then equally across reference patients."""
    from .differential import aggregate
    counts, md, _ = aggregate(ref, min_cells=1)
    cpm = counts.div(counts.sum(axis=1), axis=0) * 1e6
    cpm["cell_type"] = md.cell_type
    return cpm.groupby("cell_type", observed=True).mean().T.loc[genes]


def deconvolve_nnls(ref, target, minimum):
    """Exploratory RNA-mixture weights, NOT cell fractions or calibrated abundance.

    Patient IDs overlapping the target are excluded from reference training.
    Matching requires harmonized gene IDs; no silent symbol intersection.
    """
    overlap_patients = set(ref.obs.patient_id) & set(target.obs.patient_id)
    ref = ref[~ref.obs.patient_id.isin(overlap_patients)].copy()
    if ref.obs.patient_id.nunique() < 2:
        raise ValueError("At least two independent reference patients after overlap exclusion")
    genes = ref.var_names.intersection(target.var_names, sort=False)
    if len(genes) < minimum:
        raise ValueError("Too few shared gene IDs; harmonize the frozen genome annotation")
    signature = reference_signature(ref, genes)
    if signature.shape[1] < 2:
        raise ValueError("Reference requires at least two cell types")
    a = signature.to_numpy()
    # Downweight highly expressed genes; identical transformation on reference/target.
    scale = np.maximum(np.sqrt(a.mean(axis=1)), 1)
    a = a / scale[:, None]
    x = dense(target[:, genes].layers["counts"]).astype(float)
    all_totals = np.asarray(target.layers["counts"].sum(axis=1)).ravel()
    x = (x / all_totals[:, None] * 1e6) / scale[None, :]
    weights, residual = [], []
    for row in x:
        w, r = nnls(a, row, maxiter=100 * a.shape[1])
        if w.sum() <= 0:
            raise ValueError("Zero reference fit: out-of-reference biology or ID mismatch")
        weights.append(w / w.sum())
        residual.append(r / max(np.linalg.norm(row), 1e-12))
    return (pd.DataFrame(weights, index=target.obs_names, columns=signature.columns),
            np.asarray(residual), {"n_shared_genes": len(genes),
            "condition_number": float(np.linalg.cond(a)),
            "excluded_overlapping_patients": sorted(overlap_patients)})


def run(cfg):
    inp = Path(cfg["input_dir"]) / "spatial.h5ad"
    ref_path = Path(cfg["output_dir"]) / "singlecell/processed.h5ad"
    out = mkdir(Path(cfg["output_dir"]) / "spatial")
    a, ref = ad.read_h5ad(inp), ad.read_h5ad(ref_path)
    validate_adata(a, "spatial", cfg["mode"])
    a.X = a.layers["counts"].copy()
    x = sparse.csr_matrix(a.X)
    a.obs["total_counts"] = np.asarray(x.sum(axis=1)).ravel()
    a.obs["n_genes_by_counts"] = np.asarray((x > 0).sum(axis=1)).ravel()
    mt = a.var.gene_symbol.str.upper().str.startswith("MT-").to_numpy()
    a.obs["pct_counts_mt"] = np.asarray(x[:, mt].sum(axis=1)).ravel() / a.obs.total_counts * 100
    if not set(a.obs.in_tissue.unique()).issubset({0, 1, False, True}):
        raise ValueError("in_tissue must be binary 0/1, not arbitrary strings")
    a.obs["retained"] = ((a.obs.in_tissue == 1)
                         & (a.obs.n_genes_by_counts >= cfg["qc"]["min_spot_genes"])
                         & (a.obs.pct_counts_mt <= cfg["qc"]["max_spot_mt_pct"]))
    a.obs.to_csv(out / "spot_qc.csv")
    original_sections = set(a.obs.section_id)
    a = a[a.obs.retained].copy()
    if set(a.obs.section_id) != original_sections:
        raise ValueError("An entire section was removed by QC; review before proceeding")
    if "pathology" not in a.obs:
        if not cfg["analysis"]["allow_missing_pathology"]:
            raise ValueError("Reviewed spot pathology required for compartment analysis")
        a.obs["pathology"] = "unannotated"
    programs, coverage = scores(a, cfg["analysis"]["min_signature_coverage"])
    coverage.to_csv(out / "program_coverage.csv", index=False)
    for program in programs:
        a.obs[f"score_{program}"] = programs[program]
    method = cfg["analysis"]["deconvolution"]
    diagnostics, all_weights, morans, niches = [], [], [], []
    inputs = [inp, ref_path]
    if method == "external":
        path = Path(cfg["analysis"]["abundance_file"])
        external = pd.read_csv(path, index_col=0)
        inputs.append(path)
        if not external.index.is_unique or not external.columns.is_unique:
            raise ValueError("Duplicate external abundance IDs")
        if not set(a.obs_names).issubset(external.index):
            raise ValueError("External mapping omits retained spots")
        external = external.loc[a.obs_names]
        if not np.isfinite(external.to_numpy()).all() or (external.to_numpy() < 0).any():
            raise ValueError("Abundance must be finite and nonnegative")
        if (external.sum(axis=1) <= 0).any():
            raise ValueError("Empty abundance vector")
        external.to_csv(out / "absolute_abundance.csv")
        # Relative composition of posterior MEANS; q05 quantiles are NOT additive fractions.
        external = external.div(external.sum(axis=1), axis=0)
    elif method != "nnls":
        raise ValueError("deconvolution must be nnls or external")
    rng = np.random.default_rng(cfg["seed"])
    for section, obs in a.obs.groupby("section_id", observed=True):
        if obs.patient_id.nunique() != 1:
            raise ValueError("Each section must map to exactly one patient")
        target = a[obs.index].copy()
        w = graph(target.obsm["spatial"], cfg["analysis"]["spatial_k"])
        if method == "nnls":
            weights, residual, details = deconvolve_nnls(
                ref, target, cfg["analysis"]["min_shared_genes"])
            a.obs.loc[obs.index, "reference_residual"] = residual
            diagnostics.append({"section_id": section, **details})
        else:
            weights = external.loc[obs.index]
        all_weights.append(weights)
        for program in programs:
            val = programs.loc[obs.index, program].to_numpy()
            estimate, p = moran(val, w, cfg["analysis"]["spatial_permutations"], rng)
            morans.append({"section_id": section, "patient_id": obs.patient_id.iloc[0],
                           "program": program, "moran_I": estimate, "pvalue": p,
                           "n_spots": len(obs), "n_edges": w.nnz // 2})
        for zone in ["tumor", "stroma"]:
            mask = target.obs.pathology.eq(zone)
            if mask.sum() < 4:
                continue
            record = {"section_id": section, "patient_id": obs.patient_id.iloc[0],
                      "condition": obs.condition.iloc[0], "pathology": zone,
                      "n_spots": int(mask.sum())}
            record.update(weights.loc[mask].mean().to_dict())
            niches.append(record)
    abundance = pd.concat(all_weights).reindex(a.obs_names)
    abundance.to_csv(out / "relative_abundance.csv")
    a.obsm["relative_abundance"] = abundance
    a.uns["abundance_method"] = ("exploratory_NNLS_RNA_weights" if method == "nnls" else
                                  cfg["analysis"]["abundance_method"])
    m = pd.DataFrame(morans)
    m["qvalue"] = bh(m.pvalue)
    m.to_csv(out / "moran.csv", index=False)
    pd.DataFrame(niches, columns=["section_id", "patient_id", "condition", "pathology",
                                  "n_spots", *abundance.columns]).to_csv(
        out / "section_compartments.csv", index=False)
    write_json(out / "reference_diagnostics.json", diagnostics)
    a.write_h5ad(out / "processed.h5ad", compression="gzip")
    audit(cfg, "spatial", inputs, {"sections": a.obs.section_id.nunique(),
          "patients": a.obs.patient_id.nunique(), "mapping_method": a.uns["abundance_method"],
          "moran_null": "exchangeable spot labels within section; exploratory",
          "uncertainty": "patient bootstrap downstream; no spot-level independence assumption"})
