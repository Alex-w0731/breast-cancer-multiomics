"""Simulated counts for software verification. No patient results are generated."""

from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

from .core import MARKERS, PROGRAMS, mkdir, write_json


def generate(cfg):
    if cfg["mode"] != "demo":
        raise ValueError("Demo generation requires mode: demo")
    out = mkdir(cfg["input_dir"])
    rng = np.random.default_rng(cfg["seed"])
    genes = sorted(set(sum(MARKERS.values(), []) + sum(PROGRAMS.values(), [])))
    genes += [f"MT-DEMO{i}" for i in range(8)]
    genes += [f"DEMO{i:04d}" for i in range(800 - len(genes))]
    var = pd.DataFrame({"gene_symbol": genes}, index=pd.Index(genes, name="gene_id"))
    types = list(MARKERS)
    signature = rng.gamma(1.8, 0.8, (len(types), len(genes)))
    for i, ct in enumerate(types):
        for g in MARKERS[ct]:
            signature[i, genes.index(g)] += 18
    for g in PROGRAMS["ECM"]:
        signature[types.index("Fibroblast"), genes.index(g)] += 6
    for g in PROGRAMS["Cytotoxic"]:
        signature[types.index("T_cell"), genes.index(g)] += 6
    sc_x, sc_obs, st_x, st_obs, xy, truth, bulk, bulk_obs = [], [], [], [], [], [], [], []
    for patient in range(12):
        pid = f"SYN_P{patient:02d}"
        condition = "TNBC" if patient % 2 else "HRpos_HER2neg"
        batch = f"batch{(patient // 2) % 2}"
        sig = signature * rng.lognormal(0, 0.18, (len(types), len(genes)))
        if condition == "TNBC":
            sig[3, [genes.index(g) for g in PROGRAMS["ECM"]]] *= 1.65
        for cell in range(160):
            c = rng.choice(len(types), p=[.32, .24, .18, .18, .08])
            mu = sig[c] * rng.lognormal(0, 0.3)
            counts = rng.negative_binomial(6, 6 / (6 + mu))
            sc_x.append(counts)
            sc_obs.append({"id": f"{pid}_cell{cell:03d}", "sample_id": f"{pid}_sc",
                           "patient_id": pid, "condition": condition, "batch": batch,
                           "cell_type": types[c], "annotation_source": "synthetic_truth"})
        for spot in range(64):
            xx, yy = spot % 8, spot // 8
            tumor = (xx - 3.5) ** 2 + (yy - 3.5) ** 2 < 8
            weights = rng.dirichlet([8, 1.5, 2, 2, 1] if tumor else [1, 3, 2, 6, 1])
            mu = weights @ sig * rng.uniform(5, 12)
            st_x.append(rng.negative_binomial(15, 15 / (15 + mu)))
            sid = f"{pid}_spot{spot:03d}"
            st_obs.append({"id": sid, "sample_id": f"{pid}_st", "patient_id": pid,
                           "condition": condition, "batch": batch, "section_id": f"{pid}_s1",
                           "in_tissue": 1, "pathology": "tumor" if tumor else "stroma"})
            xy.append([xx * 100, yy * 100])
            truth.append(dict(zip(types, weights, strict=True)) | {"id": sid})
        weights = rng.dirichlet([8, 4, 3, 3, 2])
        mu = weights @ sig * 80
        bulk.append(rng.negative_binomial(20, 20 / (20 + mu)))
        bulk_obs.append({"sample_id": f"{pid}_bulk", "patient_id": pid,
                         "condition": condition, "batch": batch,
                         "age": int(rng.integers(35, 76)), "data_status": "SYNTHETIC"})
    for name, x, obs in [("singlecell", sc_x, sc_obs), ("spatial", st_x, st_obs)]:
        a = ad.AnnData(sparse.csr_matrix(np.asarray(x, dtype=np.int32)),
                       obs=pd.DataFrame(obs).set_index("id"), var=var.copy())
        a.layers["counts"] = a.X.copy()
        a.uns["data_status"] = "SYNTHETIC"
        if name == "spatial":
            a.obsm["spatial"] = np.asarray(xy, dtype=float)
        a.write_h5ad(out / f"{name}.h5ad", compression="gzip")
    pd.DataFrame(bulk, index=[r["sample_id"] for r in bulk_obs], columns=genes).to_csv(
        out / "bulk_counts.csv")
    pd.DataFrame(bulk_obs).set_index("sample_id").to_csv(out / "bulk_metadata.csv")
    var.to_csv(out / "bulk_genes.csv")
    pd.DataFrame(truth).set_index("id").to_csv(out / "spatial_truth.csv")
    write_json(out / "DATA_STATUS.json", {"status": "SYNTHETIC", "seed": cfg["seed"],
               "purpose": "software and figure demonstration only", "patients": 12})
    return Path(out)
