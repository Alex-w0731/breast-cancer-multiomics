"""Cross-modality effect tables, gene-set annotation and donor uncertainty."""

from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import hypergeom

from .core import PROGRAMS, audit, bh, mkdir, write_json


def bootstrap_mean(values, rng, n=2000):
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    if len(values) < 3:
        return np.mean(values) if len(values) else np.nan, np.nan, np.nan
    estimates = rng.choice(values, size=(n, len(values)), replace=True).mean(axis=1)
    low, high = np.quantile(estimates, [.025, .975])
    return float(values.mean()), float(low), float(high)


def ora(de):
    """Small prespecified marker sets for interpretation, not a genome-wide pathway atlas."""
    # Ambiguous symbols are excluded consistently from universe and selected hits.
    de = de.loc[~de.gene_symbol.duplicated(keep=False)].copy()
    universe = set(de.loc[de.pvalue.notna(), "gene_symbol"].dropna())
    qcol = "q_global" if "q_global" in de else "padj"
    out = []
    for direction, sign in [("up", 1), ("down", -1)]:
        hits = set(de.loc[(de[qcol] < .05) & (sign * de.log2FoldChange > 1), "gene_symbol"])
        hits &= universe
        for program, genes in PROGRAMS.items():
            gene_set = set(genes) & universe
            if not gene_set or not universe:
                continue
            k = len(gene_set & hits)
            p = hypergeom.sf(k - 1, len(universe), len(gene_set), len(hits))
            out.append({"program": program, "direction": direction, "overlap": k,
                        "set_size_tested": len(gene_set), "n_selected": len(hits),
                        "universe_size": len(universe), "pvalue": p,
                        "genes": ";".join(sorted(gene_set & hits))})
    result = pd.DataFrame(out)
    if not result.empty:
        result["qvalue"] = bh(result.pvalue)
    return result


def run(cfg):
    root = Path(cfg["output_dir"])
    out = mkdir(root / "integration")
    pb = pd.read_csv(root / "pseudobulk/differential.csv")
    bulk = pd.read_csv(root / "bulk/differential.csv")
    if bulk.gene_id.duplicated().any() or pb.duplicated(["gene_id", "cell_type"]).any():
        raise ValueError("Duplicate gene keys would inflate modality join")
    merged = pb.merge(bulk, on="gene_id", suffixes=("_sc", "_bulk"), validate="many_to_one")
    if merged.empty:
        raise ValueError("No identical gene IDs across modalities")
    merged["direction_agreement"] = (np.sign(merged.log2FoldChange_sc)
                                     == np.sign(merged.log2FoldChange_bulk))
    merged.to_csv(out / "cross_modality_effects.csv", index=False)
    enrichment = []
    for ct, sub in pb.groupby("cell_type", observed=True):
        table = ora(sub)
        table["source"] = ct
        enrichment.append(table)
    table = ora(bulk)
    table["source"] = "bulk"
    enrichment.append(table)
    enrich = pd.concat(enrichment, ignore_index=True)
    if not enrich.empty:
        enrich["q_all_program_tests"] = bh(enrich.pvalue)
    enrich.to_csv(out / "program_overrepresentation.csv", index=False)
    compartments = pd.read_csv(root / "spatial/section_compartments.csv")
    summary = []
    # Remove only this generated output, so a previous run cannot masquerade as current evidence.
    (out / "patient_tcell_compartments.csv").unlink(missing_ok=True)
    if not compartments.empty and "T_cell" in compartments:
        # Equal section weight within donor, then one tumor-stroma difference per donor.
        donor = compartments.groupby(["patient_id", "condition", "pathology"], observed=True)[
            "T_cell"].mean().unstack("pathology")
        if {"tumor", "stroma"}.issubset(donor.columns):
            donor = donor.dropna(subset=["tumor", "stroma"])
            donor["tumor_minus_stroma"] = donor.tumor - donor.stroma
            donor.to_csv(out / "patient_tcell_compartments.csv")
            rng = np.random.default_rng(cfg["seed"])
            for condition, sub in donor.reset_index().groupby("condition", observed=True):
                mean, low, high = bootstrap_mean(sub.tumor_minus_stroma, rng,
                                                cfg["analysis"]["bootstrap_replicates"])
                summary.append({"condition": condition, "n_patients": len(sub),
                                "mean_difference": mean, "ci_low": low, "ci_high": high})
    write_json(out / "spatial_patient_summary.json", summary)
    write_json(out / "interpretation_limits.json", {
        "bulk": "Subtype composition and cell-intrinsic expression are not identifiable separately",
        "cross_modality": "Direction concordance is descriptive, not a causal test or meta-analysis",
        "gene_sets": "Small author-defined marker sets; not validated clinical signatures",
        "spatial": "Relative abundance is compositional; NNLS weights reflect RNA contributions",
        "replication": "Public scRNA/spatial cohorts may overlap; TCGA is a separate cohort",
        "mode": cfg["mode"]})
    audit(cfg, "integration", [root / "pseudobulk/differential.csv", root / "bulk/differential.csv",
          root / "spatial/section_compartments.csv"], {"joined_rows": len(merged),
          "bootstrap_unit": "patient", "join": "same frozen gene_id; many cell types to one bulk gene"})
