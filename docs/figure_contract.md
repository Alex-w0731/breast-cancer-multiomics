# Figure contract / 图形规范

The delivery surface is a versioned scientific code repository with standalone
PNG (300 dpi), PDF and SVG figures. Static Matplotlib is the selected renderer.
All demo graphics carry **SYNTHETIC DATA — SOFTWARE DEMO ONLY**. Real runs are
labeled **EXPLORATORY ANALYSIS — REVIEW REQUIRED** until external review.

| Family | Question | Grain / sufficiency | Form | Interpretation boundary |
|---|---|---|---|---|
| Study design | How do three assays connect? | Proposed design; no measurements | Flow diagram | Planned work, not completed experiments |
| QC | Which observations pass prespecified rules? | Cells, spots, capture totals; >=4 observations | Histogram, scatter, bars | Descriptive; never treat cells as independent patients |
| Embedding | How do labels/batches distribute? | Retained cells | UMAP, PCA | Visualization, no patient-level significance |
| Marker matrix | Are labels consistent with marker expression? | Cell-type mean and detection fraction | Dot plot | Supporting annotation, no proof of malignancy |
| Composition | How do cell counts vary by patient? | Patient x cell type | Stacked bar / heatmap | Relative, capture-biased composition |
| Differential expression | What are adjusted condition effects? | Independent patient pseudobulks / bulk | Volcano, MA, effect matrix | Unshrunk Wald effects, explicit FDR family |
| Spatial | Where are counts, programs and estimated mixtures? | One section at a time | Equal-aspect spatial maps | No fabricated H&E; coordinate-only unless image supplied |
| Spatial autocorrelation | Are programs spatially patterned? | Section x program | Heatmap | Within-section label-permutation null |
| Patient uncertainty | What is the tumor-stroma T-cell difference? | One difference per patient, >=3 for CI | Paired points / interval | Equal section weights; patient bootstrap |
| Integration | Do condition effect directions agree? | Gene x type, shared frozen ID | Scatter / matrix | Separate cohorts, no pooling of raw assays |
| Power | What does assumed effect size imply? | Hypothetical independent patient group sizes | Curves | Planning assumptions; not achieved study power |

Use white backgrounds, dark text, quiet grids, fixed category ordering and explicit
units. Palette: blue #3569A8, gold #C69C36, orange #CE7838, olive #7B8745,
pink #B56887. Signed effects use blue/orange; continuous maps use white-to-blue.
Direct labels, marker shapes, line styles and facet titles complement color.
No decorative brand icon is applied: these are user-owned scientific figures.
Every saved figure gets a JSON manifest entry with title, input tables, data status
and formats. Numerical source tables remain beside stage outputs. A contact sheet
and Markdown gallery support final visual QA. Charts with insufficient data are
omitted with a recorded reason; no fake values or significance stars are added.

