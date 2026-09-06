"""Input contracts, counts, audit records and statistical safeguards."""

from __future__ import annotations

import hashlib
import csv
import importlib.metadata
import json
import platform
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
import yaml
from scipy import sparse
from statsmodels.stats.multitest import multipletests


PROGRAMS = {
    "ECM": ["COL1A1", "COL1A2", "COL3A1", "FN1", "POSTN", "SPARC"],
    "Cytotoxic": ["CD8A", "CD8B", "NKG7", "PRF1", "GZMB", "GNLY"],
    "IFN_response": ["ISG15", "IFIT1", "IFIT3", "MX1", "STAT1", "IRF7"],
    "Proliferation": ["MKI67", "TOP2A", "PCNA", "MCM2", "TYMS", "BIRC5"],
}
MARKERS = {
    "Epithelial": ["EPCAM", "KRT8", "KRT18", "KRT19"],
    "T_cell": ["CD3D", "CD3E", "TRAC", "CD8A"],
    "Myeloid": ["LYZ", "LST1", "TYROBP", "FCER1G"],
    "Fibroblast": ["COL1A1", "COL1A2", "DCN", "LUM"],
    "Endothelial": ["PECAM1", "VWF", "KDR", "EMCN"],
}


def config(path):
    cfg = yaml.safe_load(Path(path).read_text(encoding="utf-8"))
    if cfg["mode"] not in {"demo", "real"}:
        raise ValueError("mode must be demo or real")
    if cfg["mode"] == "real" and "demo" in Path(cfg["input_dir"]).parts:
        raise ValueError("A real run cannot use the demo input directory")
    return cfg


def mkdir(path):
    p = Path(path)
    p.mkdir(parents=True, exist_ok=True)
    return p


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(1024 * 1024), b""):
            h.update(b)
    return h.hexdigest()


def write_json(path, obj):
    def clean(value):
        if isinstance(value, dict):
            return {str(key): clean(item) for key, item in value.items()}
        if isinstance(value, (list, tuple, np.ndarray)):
            return [clean(item) for item in value]
        if isinstance(value, np.generic):
            return clean(value.item())
        if isinstance(value, float) and not np.isfinite(value):
            return None
        return value
    Path(path).write_text(json.dumps(clean(obj), indent=2, ensure_ascii=False,
                                     default=str, allow_nan=False),
                          encoding="utf-8")


def read_count_csv(path):
    # pandas renames duplicate CSV headers, which would otherwise conceal duplicate genes.
    with open(path, encoding="utf-8-sig", newline="") as f:
        header = next(csv.reader(f))
    if len(header) < 2 or len(set(header[1:])) != len(header[1:]) or any(
            not gene.strip() for gene in header[1:]):
        raise ValueError("Duplicate or blank gene IDs in count CSV header")
    return pd.read_csv(path, index_col=0)


def validate_counts(x, name="counts"):
    v = x.data if sparse.issparse(x) else np.asarray(x)
    if v.size == 0 or not np.isfinite(v).all():
        raise ValueError(f"{name}: empty or non-finite counts")
    if (v < 0).any() or not np.equal(v, np.floor(v)).all():
        raise ValueError(f"{name}: nonnegative integer raw counts required, not log/TPM")
    if np.max(v) > np.iinfo(np.int64).max:
        raise ValueError(f"{name}: integer overflow")
    if np.asarray(x.sum(axis=1)).ravel().min() <= 0:
        raise ValueError(f"{name}: zero-count observation")


def validate_adata(a, modality, mode):
    if not a.obs_names.is_unique or not a.var_names.is_unique:
        raise ValueError(f"{modality}: duplicate observation/gene IDs")
    if "counts" not in a.layers:
        raise ValueError(f"{modality}: layers['counts'] is required")
    validate_counts(a.layers["counts"], modality)
    fields = ["sample_id", "patient_id", "condition", "batch"]
    if modality == "singlecell":
        fields += ["cell_type", "annotation_source"]
    else:
        fields += ["section_id", "in_tissue"]
        if "spatial" not in a.obsm or a.obsm["spatial"].shape != (a.n_obs, 2):
            raise ValueError("spatial: finite full-resolution x/y coordinates required")
        if not np.isfinite(a.obsm["spatial"]).all():
            raise ValueError("spatial: non-finite coordinates")
    for col in fields:
        if col not in a.obs or a.obs[col].isna().any():
            raise ValueError(f"{modality}: missing metadata: {col}")
        if (a.obs[col].astype(str).str.strip() == "").any():
            raise ValueError(f"{modality}: blank metadata: {col}")
    for col in ["patient_id", "condition", "batch"]:
        if a.obs.groupby("sample_id", observed=True)[col].nunique().max() != 1:
            raise ValueError(f"{modality}: each sample must map to exactly one {col}")
    if a.obs.groupby("patient_id", observed=True).condition.nunique().max() != 1:
        raise ValueError("This between-patient design requires one condition per patient")
    if "gene_symbol" not in a.var:
        raise ValueError("var['gene_symbol'] must be supplied from the frozen annotation")
    if a.var.gene_symbol.isna().any():
        raise ValueError("Missing gene symbols: curate the reference mapping first")
    if mode == "real" and a.uns.get("data_status") == "SYNTHETIC":
        raise ValueError("Synthetic data cannot enter a real run")


def dense(x):
    return x.toarray() if sparse.issparse(x) else np.asarray(x)


def log_cpm(x):
    x = sparse.csr_matrix(x, dtype=np.float64)
    totals = np.asarray(x.sum(axis=1)).ravel()
    if (totals <= 0).any():
        raise ValueError("Cannot normalize empty libraries")
    out = sparse.diags(1e6 / totals) @ x
    out.data = np.log1p(out.data)
    return out.tocsr()


def bh(p):
    p = np.asarray(p, dtype=float)
    valid = np.isfinite(p)
    q = np.full(p.shape, np.nan)
    if valid.any():
        q[valid] = multipletests(p[valid], method="fdr_bh")[1]
    return q


def scores(a, coverage=0.6):
    """Descriptive mean log-CPM, no cross-platform absolute comparisons.

    One-to-many gene symbols are excluded to avoid double weighting.
    Raw Ensembl IDs remain the feature key for DE and modality joins.
    """
    x = log_cpm(a.layers["counts"])
    symbols = a.var.gene_symbol.astype(str)
    unique = ~symbols.duplicated(keep=False)
    out, audit = {}, []
    for name, genes in PROGRAMS.items():
        mask = symbols.isin(genes) & unique
        found = symbols[mask].tolist()
        fraction = len(found) / len(genes)
        audit.append({"program": name, "coverage": fraction, "genes": ";".join(found)})
        if fraction < coverage:
            raise ValueError(f"{name}: marker coverage {fraction:.0%} < {coverage:.0%}")
        out[name] = np.asarray(x[:, mask.to_numpy()].mean(axis=1)).ravel()
    return pd.DataFrame(out, index=a.obs_names), pd.DataFrame(audit)


def audit(cfg, stage, inputs, extra=None):
    out = mkdir(Path(cfg["output_dir"]) / "audit")
    packages = {}
    for package in ["bcflow", "numpy", "pandas", "scanpy", "anndata", "pydeseq2", "scipy"]:
        try:
            packages[package] = importlib.metadata.version(package)
        except importlib.metadata.PackageNotFoundError:
            packages[package] = "not-installed"
    record = {"stage": stage, "utc": datetime.now(timezone.utc).isoformat(),
              "mode": cfg["mode"], "seed": cfg["seed"], "config": cfg,
              "python": platform.python_version(), "platform": platform.platform(),
              "packages": packages,
              "inputs": {str(p): sha256(p) for p in inputs}, "details": extra or {}}
    record["source_sha256"] = {p.name: sha256(p) for p in Path(__file__).parent.glob("*.py")}
    write_json(out / f"{stage}.json", record)
