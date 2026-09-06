"""Publication-oriented exports with explicit data status and traceable source tables."""

from pathlib import Path
import json

import anndata as ad
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch
from PIL import Image, ImageDraw
from sklearn.decomposition import PCA
from statsmodels.stats.power import TTestIndPower

from .core import MARKERS, PROGRAMS, dense, mkdir, write_json


COLORS = ["#3569A8", "#C69C36", "#CE7838", "#7B8745", "#B56887"]
CMAP = LinearSegmentedColormap.from_list("bc_blue", ["#F4F7FA", "#3569A8"])
DIVERGING = LinearSegmentedColormap.from_list("bc_signed", [COLORS[0], "#FAFAFA", COLORS[2]])


class Gallery:
    def __init__(self, cfg):
        self.cfg = cfg
        self.out = mkdir(Path(cfg["output_dir"]) / "figures")
        self.entries = []
        self.omissions = []
        self.status = ("SYNTHETIC DATA - SOFTWARE DEMO ONLY" if cfg["mode"] == "demo" else
                       "EXPLORATORY ANALYSIS - REVIEW REQUIRED")
        plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10,
                             "axes.spines.top": False, "axes.spines.right": False,
                             "axes.titleweight": "semibold", "pdf.fonttype": 42,
                             "svg.fonttype": "none", "savefig.facecolor": "white"})

    def save(self, fig, slug, title, subtitle, source):
        fig.suptitle(title, x=.06, y=.99, ha="left", fontsize=14, fontweight="bold")
        fig.text(.06, .94, subtitle, ha="left", va="top", fontsize=9, color="#454545")
        fig.text(.06, .025, self.status, fontsize=8, color="#77552D", fontweight="bold")
        fig.tight_layout(rect=(.025, .065, .97, .86))
        for ext in ["png", "pdf", "svg"]:
            fig.savefig(self.out / f"{slug}.{ext}", dpi=300)
        preview_dir = mkdir(self.out / "previews")
        fig.savefig(preview_dir / f"{slug}.png", dpi=120)
        self.entries.append({"id": slug, "title": title, "subtitle": subtitle,
                             "status": self.status, "sources": source,
                             "files": [f"{slug}.{e}" for e in ["png", "pdf", "svg"]]})
        plt.close(fig)

    def finish(self):
        write_json(self.out / "figure_manifest.json", self.entries)
        write_json(self.out / "omitted_figures.json", self.omissions)
        md = ["# Figure gallery", "", f"**{self.status}**", "",
              "Each PNG has matching editable SVG and vector PDF exports.", ""]
        thumb_w, thumb_h = 480, 350
        sheet = Image.new("RGB", (4 * thumb_w, ((len(self.entries) + 3) // 4) * thumb_h), "white")
        draw = ImageDraw.Draw(sheet)
        for i, e in enumerate(self.entries):
            md += [f"## {e['id']} — {e['title']}", "", e["subtitle"], "",
                   f"![{e['title']}]({e['id']}.png)", ""]
            with Image.open(self.out / f"{e['id']}.png") as im:
                im.thumbnail((thumb_w - 12, thumb_h - 25))
                x, y = (i % 4) * thumb_w, (i // 4) * thumb_h
                sheet.paste(im, (x, y + 20))
                draw.text((x + 10, y + 3), e["id"], fill="#333333")
        sheet.save(self.out / "contact_sheet.jpg", quality=90)
        (self.out / "GALLERY.md").write_text("\n".join(md), encoding="utf-8")


def embedding(g, a, key, slug, title):
    fig, ax = plt.subplots(figsize=(8.4, 6.6))
    values = a.obs[key].astype(str)
    cats = sorted(values.unique())
    markers = ["o", "s", "^", "D", "v", "P", "X", "<", ">", "h", "p", "*"]
    for i, cat in enumerate(cats):
        mask = values.eq(cat)
        coords = a.obsm["X_umap"][mask]
        ax.scatter(coords[:, 0], coords[:, 1], s=7, alpha=.6, linewidth=0,
                   color=COLORS[i % 5], marker=markers[i % len(markers)], label=cat,
                   rasterized=True)
    ax.set(xlabel="UMAP 1", ylabel="UMAP 2")
    ax.legend(loc="center left", bbox_to_anchor=(1.01, .5), frameon=False, fontsize=8,
              markerscale=2)
    g.save(fig, slug, title, f"{a.n_obs:,} retained cells; labels are descriptive",
           ["singlecell/processed.h5ad"])


def volcano(g, table, slug, title, qcol, source):
    fig, ax = plt.subplots(figsize=(8, 6.2))
    valid = table.dropna(subset=["log2FoldChange", qcol])
    x, q = valid.log2FoldChange, valid[qcol]
    for name, mask, col in [
        ("Other tested", (q >= .05) | (x.abs() <= 1), "#C3C7CA"),
        ("Lower, q < .05", (q < .05) & (x < -1), COLORS[0]),
        ("Higher, q < .05", (q < .05) & (x > 1), COLORS[2]),
    ]:
        ax.scatter(x[mask], -np.log10(q[mask].clip(lower=1e-300)), s=13, color=col,
                   alpha=.75, label=name, rasterized=True)
    ax.axhline(-np.log10(.05), color="#555555", ls="--", lw=.8)
    ax.axvline(0, color="#888888", lw=.6)
    ax.set(xlabel="Unshrunk log2 fold change: TNBC / HR+ HER2-", ylabel=f"-log10({qcol})")
    ax.legend(frameon=False, fontsize=8)
    g.save(fig, slug, title, f"{len(valid):,} genes with finite adjusted p; {qcol} controls testing family", source)


def run(cfg):
    root = Path(cfg["output_dir"])
    g = Gallery(cfg)
    a = ad.read_h5ad(root / "singlecell/processed.h5ad")
    st = ad.read_h5ad(root / "spatial/processed.h5ad")
    qc = pd.read_csv(root / "singlecell/cell_qc.csv", index_col=0)
    bulk = pd.read_csv(root / "bulk/differential.csv")
    pb = pd.read_csv(root / "pseudobulk/differential.csv")
    # F01: actual workflow structure, not invented quantitative evidence.
    fig, ax = plt.subplots(figsize=(11, 5.8))
    ax.set(xlim=(0, 12), ylim=(0, 5))
    ax.axis("off")
    boxes = [(0, 2, "Patient + pathology\nprospective paired sampling"),
             (3.5, 3.7, "scRNA-seq\nQC / labels / pseudobulk"),
             (3.5, 2, "Spatial RNA\nQC / mapping / niches"),
             (3.5, .3, "Bulk RNA-seq\ncounts / adjusted DE"),
             (7.7, 2, "Cross-modal evidence\npatient-level validation")]
    for i, (x, y, label) in enumerate(boxes):
        ax.add_patch(FancyBboxPatch((x + .1, y), 3, 1, boxstyle="round,pad=.08",
                                   facecolor="#F3F6F9", edgecolor=COLORS[i % 5]))
        ax.text(x + 1.6, y + .5, label, ha="center", va="center", fontsize=10)
    for y in [4.2, 2.5, .8]:
        ax.add_patch(FancyArrowPatch((3.2, 2.5), (3.5, y), arrowstyle="->", mutation_scale=12))
        ax.add_patch(FancyArrowPatch((6.7, y), (7.7, 2.5), arrowstyle="->", mutation_scale=12))
    g.save(fig, "F01_workflow", "Three-modality analysis design",
           "Proposed study structure; donor identity is preserved across every assay", ["docs/experimental_design_zh.md"])
    # QC distributions and library/complexity relationship.
    for number, col, label in [(2, "total_counts", "Total UMI counts / cell"),
                                (3, "n_genes_by_counts", "Detected genes / cell"),
                                (4, "pct_counts_mt", "Mitochondrial counts (%)")]:
        fig, ax = plt.subplots(figsize=(8, 5.7))
        ax.hist(qc[col], bins=45, color=COLORS[0], edgecolor="white", linewidth=.3)
        ax.set(xlabel=label, ylabel="Number of cells")
        g.save(fig, f"F{number:02d}_qc_{col}", f"Single-cell QC: {label}",
               f"{len(qc):,} input cells before filtering; counts are untransformed", ["singlecell/cell_qc.csv"])
    fig, ax = plt.subplots(figsize=(8, 6))
    points = ax.scatter(qc.total_counts, qc.n_genes_by_counts, c=qc.pct_counts_mt,
                        cmap=CMAP, s=7, rasterized=True)
    fig.colorbar(points, ax=ax, label="Mitochondrial counts (%)")
    ax.set(xlabel="Total UMI counts / cell", ylabel="Detected genes / cell")
    g.save(fig, "F05_complexity", "Library complexity and mitochondrial fraction",
           f"{len(qc):,} input cells; each point is one cell", ["singlecell/cell_qc.csv"])
    fig, ax = plt.subplots(figsize=(9, 5.8))
    retention = qc.groupby("sample_id").retained.agg(["size", "sum"])
    ax.bar(np.arange(len(retention)), retention["size"], color="#DADDE0", label="Input")
    ax.bar(np.arange(len(retention)), retention["sum"], color=COLORS[0], label="Retained")
    ax.set(xticks=np.arange(len(retention)), xticklabels=retention.index, ylabel="Cells")
    ax.tick_params(axis="x", rotation=60)
    ax.legend(frameon=False)
    g.save(fig, "F06_retention", "Cell retention by capture", "Input and retained cells; capture is a technical unit", ["singlecell/cell_qc.csv"])
    for i, key in enumerate(["cell_type", "condition", "batch", "leiden"], 7):
        embedding(g, a, key, f"F{i:02d}_umap_{key}", f"Single-cell UMAP by {key}")
    fig, ax = plt.subplots(figsize=(8, 5.7))
    variance = a.uns["pca"]["variance_ratio"]
    ax.plot(np.arange(1, len(variance) + 1), variance * 100, "o-", color=COLORS[0], ms=4)
    ax.set(xlabel="Principal component", ylabel="Explained variance (%)")
    g.save(fig, "F11_pca_variance", "PCA variance explained", "HVG-based uncorrected embedding; inspect batch effects separately", ["singlecell/processed.h5ad"])
    fig, ax = plt.subplots(figsize=(11, 6))
    marker_genes = [x for x in sum(MARKERS.values(), []) if x in set(a.var.gene_symbol)]
    cell_types = sorted(a.obs.cell_type.astype(str).unique())
    dots = []
    for i, ct in enumerate(cell_types):
        mask = a.obs.cell_type.eq(ct)
        for j, gene in enumerate(marker_genes):
            subset = a[mask, a.var.gene_symbol.eq(gene)].X
            mean = float(np.mean(dense(subset)))
            fraction = float(np.mean(dense(subset) > 0))
            dots.append((j, i, fraction, mean))
    dots = np.asarray(dots)
    point = ax.scatter(dots[:, 0], dots[:, 1], s=90 * dots[:, 2], c=dots[:, 3],
                       cmap=CMAP, vmin=0, vmax=max(float(dots[:, 3].max()), 1e-6),
                       edgecolors="#777777", linewidths=.3)
    handles = [ax.scatter([], [], s=90 * fraction, color="#8DA9C8", edgecolors="#777777",
                          linewidths=.3, label=f"{fraction:.0%}")
               for fraction in [.25, .5, .75, 1.0]]
    ax.legend(handles=handles, title="Fraction with detected expression", ncol=4,
              loc="lower left", bbox_to_anchor=(0, 1.01), frameon=False,
              fontsize=8, title_fontsize=8)
    fig.colorbar(point, ax=ax, label="Mean log1p(normalized counts)")
    ax.set(xticks=range(len(marker_genes)), xticklabels=marker_genes,
           yticks=range(len(cell_types)), yticklabels=cell_types)
    ax.tick_params(axis="x", rotation=75)
    g.save(fig, "F12_marker_dotplot", "Annotation marker expression",
           "Dot area scales with detection fraction; marker support does not establish malignancy", ["singlecell/processed.h5ad"])
    composition = pd.read_csv(root / "singlecell/patient_composition.csv", index_col=0)
    fig, ax = plt.subplots(figsize=(10, 6))
    composition.plot.bar(stacked=True, ax=ax, color=COLORS, width=.8)
    ax.set(ylabel="Fraction of retained cells", xlabel="Patient", ylim=(0, 1))
    ax.legend(loc="center left", bbox_to_anchor=(1, .5), frameon=False)
    g.save(fig, "F13_composition", "Cell composition by patient", "Denominator: all retained cells from each patient", ["singlecell/patient_composition.csv"])
    fig, ax = plt.subplots(figsize=(8.5, 6.5))
    sns.heatmap(composition.T, cmap=CMAP, vmin=0, vmax=1, ax=ax, cbar_kws={"label": "Cell fraction"})
    g.save(fig, "F14_composition_matrix", "Patient-by-cell-type composition", "Relative abundance reflects capture and dissociation biases", ["singlecell/patient_composition.csv"])
    # Patient-level program distributions: each point is one patient, one cell type.
    programs = pd.read_csv(root / "singlecell/patient_programs.csv")
    for number, program in enumerate(PROGRAMS, 15):
        fig, ax = plt.subplots(figsize=(10, 6))
        sns.boxplot(data=programs, x="cell_type", y=f"score_{program}", hue="condition", ax=ax,
                    palette=COLORS[:2], fliersize=0)
        sns.stripplot(data=programs, x="cell_type", y=f"score_{program}", hue="condition", ax=ax,
                      palette=COLORS[:2], dodge=True, jitter=False, size=4, edgecolor="#333333",
                      linewidth=.4, legend=False)
        ax.set(xlabel="Reviewed cell type", ylabel="Mean marker log1p(CPM)")
        ax.legend(frameon=False, fontsize=8)
        g.save(fig, f"F{number:02d}_program_{program}", f"Patient-level {program} marker score",
               "One point per patient and type; scores are descriptive within this assay", ["singlecell/patient_programs.csv"])
    for number, (ct, sub) in enumerate(pb.groupby("cell_type", observed=True), 19):
        volcano(g, sub, f"F{number:02d}_pseudobulk_{ct}", f"Patient-pseudobulk DE: {ct}",
                "q_global", ["pseudobulk/differential.csv", "pseudobulk/metadata.csv"])
    # Fixed B/S/I prefixes keep IDs unique if real data have more than five cell types.
    volcano(g, bulk, "B01_bulk_volcano", "Bulk RNA-seq differential expression", "padj", ["bulk/differential.csv"])
    fig, ax = plt.subplots(figsize=(8, 6))
    valid = bulk.dropna(subset=["baseMean", "log2FoldChange"])
    ax.scatter(np.log10(valid.baseMean + 1), valid.log2FoldChange, s=10,
               c=np.where(valid.padj < .05, COLORS[2], "#C7CBCF"), rasterized=True)
    ax.axhline(0, color="#555555", lw=.8)
    ax.set(xlabel="log10(mean normalized count + 1)", ylabel="Unshrunk log2 fold change")
    g.save(fig, "B02_bulk_MA", "Bulk RNA-seq mean-expression versus effect", "Orange: BH padj < .05; effect estimates are unshrunk", ["bulk/differential.csv"])
    normalized = pd.read_csv(root / "bulk/normalized_counts.csv", index_col=0)
    bmd = pd.read_csv(root / "bulk/metadata.csv", index_col=0).loc[normalized.index]
    lognorm = np.log1p(normalized)
    bpc = PCA(n_components=2, random_state=cfg["seed"]).fit_transform(lognorm)
    fig, ax = plt.subplots(figsize=(8, 6))
    for j, cond in enumerate(sorted(bmd.condition.unique())):
        mask = bmd.condition.eq(cond)
        ax.scatter(bpc[mask, 0], bpc[mask, 1], c=COLORS[j], marker=["o", "s"][j], label=cond)
    ax.legend(frameon=False)
    ax.set(xlabel="PC 1", ylabel="PC 2")
    g.save(fig, "B03_bulk_PCA", "Bulk sample PCA", f"{len(bmd)} independent patients; log1p size-factor-normalized counts", ["bulk/normalized_counts.csv", "bulk/metadata.csv"])
    fig, ax = plt.subplots(figsize=(8.5, 7))
    sns.heatmap(lognorm.T.corr(method="spearman"), ax=ax, cmap=CMAP, vmin=0, vmax=1)
    g.save(fig, "B04_bulk_correlation", "Bulk sample expression correlation", "Spearman correlation across expressed genes; descriptive QC", ["bulk/normalized_counts.csv"])
    fig, ax = plt.subplots(figsize=(10, 6.5))
    variable = lognorm.var().nlargest(25).index
    z = (lognorm[variable] - lognorm[variable].mean()) / lognorm[variable].std().replace(0, 1)
    sns.heatmap(z.T, ax=ax, cmap=DIVERGING, center=0, xticklabels=True, yticklabels=True,
                cbar_kws={"label": "Within-gene z score"})
    g.save(fig, "B05_bulk_variable_genes", "Most variable bulk expression features", "25 genes chosen by variance, not by favorable differential-expression p values", ["bulk/normalized_counts.csv"])
    # One explicit section selected by stable ID. All sections remain in source objects/tables.
    section = sorted(st.obs.section_id.astype(str).unique())[0]
    s = st[st.obs.section_id == section].copy()
    maps = {"UMI_counts": s.obs.total_counts, "Detected_genes": s.obs.n_genes_by_counts,
            "Mitochondrial_pct": s.obs.pct_counts_mt}
    maps.update({p: s.obs[f"score_{p}"] for p in PROGRAMS})
    abundance = pd.read_csv(root / "spatial/relative_abundance.csv", index_col=0).loc[s.obs_names]
    maps.update({f"Mixture_{ct}": abundance[ct] for ct in abundance})
    for number, (name, values) in enumerate(maps.items(), 1):
        fig, ax = plt.subplots(figsize=(8, 7))
        coords = s.obsm["spatial"]
        pts = ax.scatter(coords[:, 0], coords[:, 1], c=values, cmap=CMAP, s=36,
                         edgecolor="#666666", linewidth=.2, rasterized=True)
        ax.set_aspect("equal")
        ax.invert_yaxis()
        ax.set(xlabel="Full-resolution image x (pixels)", ylabel="Full-resolution image y (pixels)")
        fig.colorbar(pts, ax=ax, label=name.replace("_", " "))
        g.save(fig, f"S{number:02d}_{name}", f"Spatial map: {name.replace('_', ' ')}",
               f"Section {section}; {s.n_obs} retained spots; coordinate-only map", ["spatial/processed.h5ad", "spatial/relative_abundance.csv"])
    fig, ax = plt.subplots(figsize=(8.5, 6.8))
    mor = pd.read_csv(root / "spatial/moran.csv")
    sns.heatmap(mor.pivot(index="section_id", columns="program", values="moran_I"), ax=ax,
                cmap=DIVERGING, center=0, annot=True, fmt=".2f", cbar_kws={"label": "Moran's I"})
    g.save(fig, "S20_moran", "Section-level spatial autocorrelation", "Each graph is constructed within one section; inferential details are in moran.csv", ["spatial/moran.csv"])
    niche_path = root / "integration/patient_tcell_compartments.csv"
    if niche_path.exists():
        niche = pd.read_csv(niche_path)
        fig, ax = plt.subplots(figsize=(8, 6))
        for _, row in niche.iterrows():
            ax.plot([0, 1], [row.stroma, row.tumor], color="#858A91", alpha=.65,
                    marker="o", lw=1, markerfacecolor="white")
        ax.set(xticks=[0, 1], xticklabels=["Stroma", "Tumor"], ylabel="Relative T-cell mixture weight")
        g.save(fig, "S21_patient_compartments", "T-cell mixture by tissue compartment",
               f"{len(niche)} paired patients; sections averaged equally within each patient", ["integration/patient_tcell_compartments.csv"])
    else:
        g.omissions.append({"id": "S21", "reason": "No evaluable paired pathology compartments / T_cell column"})
    intervals = pd.DataFrame(json.loads((root / "integration/spatial_patient_summary.json").read_text()))
    if not intervals.empty and intervals[["ci_low", "ci_high"]].notna().all().all():
        fig, ax = plt.subplots(figsize=(9.5, max(4.0, 2.8 + .5 * len(intervals))))
        for i, row in intervals.iterrows():
            ax.hlines(i, row.ci_low, row.ci_high, color=COLORS[0], lw=1.5)
            ax.plot([row.ci_low, row.ci_high], [i, i], "|", color=COLORS[0], ms=8)
            ax.plot(row.mean_difference, i, "o", color=COLORS[0])
        ax.axvline(0, color="#555555", ls="--", lw=.8)
        ax.set(yticks=range(len(intervals)), yticklabels=[f"{r.condition} (n={r.n_patients})"
               for r in intervals.itertuples()], xlabel="Mean tumor minus stroma T-cell mixture weight",
               ylim=(-.5, len(intervals) - .5))
        g.save(fig, "S22_patient_uncertainty", "Patient-level compartment differences",
               "95% percentile bootstrap intervals; resampling unit is patient, not spot",
               ["integration/spatial_patient_summary.json", "integration/patient_tcell_compartments.csv"])
    merged = pd.read_csv(root / "integration/cross_modality_effects.csv")
    ct = "Fibroblast" if "Fibroblast" in merged.cell_type.unique() else merged.cell_type.iloc[0]
    sub = merged[merged.cell_type == ct].dropna(subset=["log2FoldChange_sc", "log2FoldChange_bulk"])
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.scatter(sub.log2FoldChange_sc, sub.log2FoldChange_bulk, s=10, alpha=.55,
               color=COLORS[0], rasterized=True)
    ax.axhline(0, color="#777777", lw=.8)
    ax.axvline(0, color="#777777", lw=.8)
    ax.set(xlabel=f"{ct} pseudobulk log2FC", ylabel="Bulk log2FC")
    g.save(fig, "I01_effect_concordance", "Condition effects across modalities",
           f"{len(sub)} shared gene IDs; agreement does not separate composition from regulation", ["integration/cross_modality_effects.csv"])
    fig, ax = plt.subplots(figsize=(9, 6))
    ecm = merged[merged.gene_symbol_sc.isin(PROGRAMS["ECM"])].pivot(
        index="gene_symbol_sc", columns="cell_type", values="log2FoldChange_sc")
    sns.heatmap(ecm, ax=ax, cmap=DIVERGING, center=0, annot=True, fmt=".2f",
                cbar_kws={"label": "Unshrunk log2FC"})
    g.save(fig, "I02_ECM_effects", "ECM marker effects across cell types", "Prespecified marker set; gene-level statistics remain in the full source table", ["pseudobulk/differential.csv"])
    fig, ax = plt.subplots(figsize=(8, 6))
    power = TTestIndPower()
    group_n = np.arange(6, 81, 2)
    for i, effect in enumerate([.4, .6, .8]):
        ax.plot(group_n, power.power(effect_size=effect, nobs1=group_n, alpha=.05,
                                     ratio=1, alternative="two-sided"),
                color=COLORS[i], ls=["-", "--", ":"][i], label=f"Cohen d = {effect}")
    ax.axhline(.8, color="#555555", ls="--", lw=.8)
    ax.set(xlabel="Independent patients per group", ylabel="Theoretical power", ylim=(0, 1))
    ax.legend(frameon=False)
    g.save(fig, "D01_power_sensitivity", "Patient sample-size sensitivity",
           "Hypothetical two-sample t test; alpha=.05; no dropout, multiplicity or batch penalty", ["docs/experimental_design_zh.md"])
    g.finish()
    return len(g.entries)
