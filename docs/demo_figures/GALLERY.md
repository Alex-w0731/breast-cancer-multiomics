# Figure gallery

**SYNTHETIC DATA - SOFTWARE DEMO ONLY**

This gallery contains 46 figures. Repository PNGs are 120 dpi browsing previews.
The downloadable local delivery and successful GitHub Actions artifacts contain
all 138 full exports: 300 dpi PNG, editable SVG and vector PDF for each figure.
Rerunning the pipeline regenerates these exports in `results/demo/figures/`.

## F01_workflow — Three-modality analysis design

Proposed study structure; donor identity is preserved across every assay

![Three-modality analysis design](F01_workflow.png)

## F02_qc_total_counts — Single-cell QC: Total UMI counts / cell

1,920 input cells before filtering; counts are untransformed

![Single-cell QC: Total UMI counts / cell](F02_qc_total_counts.png)

## F03_qc_n_genes_by_counts — Single-cell QC: Detected genes / cell

1,920 input cells before filtering; counts are untransformed

![Single-cell QC: Detected genes / cell](F03_qc_n_genes_by_counts.png)

## F04_qc_pct_counts_mt — Single-cell QC: Mitochondrial counts (%)

1,920 input cells before filtering; counts are untransformed

![Single-cell QC: Mitochondrial counts (%)](F04_qc_pct_counts_mt.png)

## F05_complexity — Library complexity and mitochondrial fraction

1,920 input cells; each point is one cell

![Library complexity and mitochondrial fraction](F05_complexity.png)

## F06_retention — Cell retention by capture

Input and retained cells; capture is a technical unit

![Cell retention by capture](F06_retention.png)

## F07_umap_cell_type — Single-cell UMAP by cell_type

1,920 retained cells; labels are descriptive

![Single-cell UMAP by cell_type](F07_umap_cell_type.png)

## F08_umap_condition — Single-cell UMAP by condition

1,920 retained cells; labels are descriptive

![Single-cell UMAP by condition](F08_umap_condition.png)

## F09_umap_batch — Single-cell UMAP by batch

1,920 retained cells; labels are descriptive

![Single-cell UMAP by batch](F09_umap_batch.png)

## F10_umap_leiden — Single-cell UMAP by leiden

1,920 retained cells; labels are descriptive

![Single-cell UMAP by leiden](F10_umap_leiden.png)

## F11_pca_variance — PCA variance explained

HVG-based uncorrected embedding; inspect batch effects separately

![PCA variance explained](F11_pca_variance.png)

## F12_marker_dotplot — Annotation marker expression

Dot area scales with detection fraction; marker support does not establish malignancy

![Annotation marker expression](F12_marker_dotplot.png)

## F13_composition — Cell composition by patient

Denominator: all retained cells from each patient

![Cell composition by patient](F13_composition.png)

## F14_composition_matrix — Patient-by-cell-type composition

Relative abundance reflects capture and dissociation biases

![Patient-by-cell-type composition](F14_composition_matrix.png)

## F15_program_ECM — Patient-level ECM marker score

One point per patient and type; scores are descriptive within this assay

![Patient-level ECM marker score](F15_program_ECM.png)

## F16_program_Cytotoxic — Patient-level Cytotoxic marker score

One point per patient and type; scores are descriptive within this assay

![Patient-level Cytotoxic marker score](F16_program_Cytotoxic.png)

## F17_program_IFN_response — Patient-level IFN_response marker score

One point per patient and type; scores are descriptive within this assay

![Patient-level IFN_response marker score](F17_program_IFN_response.png)

## F18_program_Proliferation — Patient-level Proliferation marker score

One point per patient and type; scores are descriptive within this assay

![Patient-level Proliferation marker score](F18_program_Proliferation.png)

## F19_pseudobulk_Endothelial — Patient-pseudobulk DE: Endothelial

664 genes with finite adjusted p; q_global controls testing family

![Patient-pseudobulk DE: Endothelial](F19_pseudobulk_Endothelial.png)

## F20_pseudobulk_Epithelial — Patient-pseudobulk DE: Epithelial

785 genes with finite adjusted p; q_global controls testing family

![Patient-pseudobulk DE: Epithelial](F20_pseudobulk_Epithelial.png)

## F21_pseudobulk_Fibroblast — Patient-pseudobulk DE: Fibroblast

754 genes with finite adjusted p; q_global controls testing family

![Patient-pseudobulk DE: Fibroblast](F21_pseudobulk_Fibroblast.png)

## F22_pseudobulk_Myeloid — Patient-pseudobulk DE: Myeloid

759 genes with finite adjusted p; q_global controls testing family

![Patient-pseudobulk DE: Myeloid](F22_pseudobulk_Myeloid.png)

## F23_pseudobulk_T_cell — Patient-pseudobulk DE: T_cell

770 genes with finite adjusted p; q_global controls testing family

![Patient-pseudobulk DE: T_cell](F23_pseudobulk_T_cell.png)

## B01_bulk_volcano — Bulk RNA-seq differential expression

799 genes with finite adjusted p; padj controls testing family

![Bulk RNA-seq differential expression](B01_bulk_volcano.png)

## B02_bulk_MA — Bulk RNA-seq mean-expression versus effect

Orange: BH padj < .05; effect estimates are unshrunk

![Bulk RNA-seq mean-expression versus effect](B02_bulk_MA.png)

## B03_bulk_PCA — Bulk sample PCA

12 independent patients; log1p size-factor-normalized counts

![Bulk sample PCA](B03_bulk_PCA.png)

## B04_bulk_correlation — Bulk sample expression correlation

Spearman correlation across expressed genes; descriptive QC

![Bulk sample expression correlation](B04_bulk_correlation.png)

## B05_bulk_variable_genes — Most variable bulk expression features

25 genes chosen by variance, not by favorable differential-expression p values

![Most variable bulk expression features](B05_bulk_variable_genes.png)

## S01_UMI_counts — Spatial map: UMI counts

Section SYN_P00_s1; 64 retained spots; coordinate-only map

![Spatial map: UMI counts](S01_UMI_counts.png)

## S02_Detected_genes — Spatial map: Detected genes

Section SYN_P00_s1; 64 retained spots; coordinate-only map

![Spatial map: Detected genes](S02_Detected_genes.png)

## S03_Mitochondrial_pct — Spatial map: Mitochondrial pct

Section SYN_P00_s1; 64 retained spots; coordinate-only map

![Spatial map: Mitochondrial pct](S03_Mitochondrial_pct.png)

## S04_ECM — Spatial map: ECM

Section SYN_P00_s1; 64 retained spots; coordinate-only map

![Spatial map: ECM](S04_ECM.png)

## S05_Cytotoxic — Spatial map: Cytotoxic

Section SYN_P00_s1; 64 retained spots; coordinate-only map

![Spatial map: Cytotoxic](S05_Cytotoxic.png)

## S06_IFN_response — Spatial map: IFN response

Section SYN_P00_s1; 64 retained spots; coordinate-only map

![Spatial map: IFN response](S06_IFN_response.png)

## S07_Proliferation — Spatial map: Proliferation

Section SYN_P00_s1; 64 retained spots; coordinate-only map

![Spatial map: Proliferation](S07_Proliferation.png)

## S08_Mixture_Endothelial — Spatial map: Mixture Endothelial

Section SYN_P00_s1; 64 retained spots; coordinate-only map

![Spatial map: Mixture Endothelial](S08_Mixture_Endothelial.png)

## S09_Mixture_Epithelial — Spatial map: Mixture Epithelial

Section SYN_P00_s1; 64 retained spots; coordinate-only map

![Spatial map: Mixture Epithelial](S09_Mixture_Epithelial.png)

## S10_Mixture_Fibroblast — Spatial map: Mixture Fibroblast

Section SYN_P00_s1; 64 retained spots; coordinate-only map

![Spatial map: Mixture Fibroblast](S10_Mixture_Fibroblast.png)

## S11_Mixture_Myeloid — Spatial map: Mixture Myeloid

Section SYN_P00_s1; 64 retained spots; coordinate-only map

![Spatial map: Mixture Myeloid](S11_Mixture_Myeloid.png)

## S12_Mixture_T_cell — Spatial map: Mixture T cell

Section SYN_P00_s1; 64 retained spots; coordinate-only map

![Spatial map: Mixture T cell](S12_Mixture_T_cell.png)

## S20_moran — Section-level spatial autocorrelation

Each graph is constructed within one section; inferential details are in moran.csv

![Section-level spatial autocorrelation](S20_moran.png)

## S21_patient_compartments — T-cell mixture by tissue compartment

12 paired patients; sections averaged equally within each patient

![T-cell mixture by tissue compartment](S21_patient_compartments.png)

## S22_patient_uncertainty — Patient-level compartment differences

95% percentile bootstrap intervals; resampling unit is patient, not spot

![Patient-level compartment differences](S22_patient_uncertainty.png)

## I01_effect_concordance — Condition effects across modalities

754 shared gene IDs; agreement does not separate composition from regulation

![Condition effects across modalities](I01_effect_concordance.png)

## I02_ECM_effects — ECM marker effects across cell types

Prespecified marker set; gene-level statistics remain in the full source table

![ECM marker effects across cell types](I02_ECM_effects.png)

## D01_power_sensitivity — Patient sample-size sensitivity

Hypothetical two-sample t test; alpha=.05; no dropout, multiplicity or batch penalty

![Patient sample-size sensitivity](D01_power_sensitivity.png)
