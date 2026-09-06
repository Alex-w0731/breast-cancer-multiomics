# Statistical analysis plan (SAP), v0.1

**Status: prospective specification and implemented exploratory workflow. No real-cohort result is asserted.**

## Population and estimands

Breast primary invasive carcinomas, prespecified untreated cohort. Compare clinically
verified TNBC against HR-positive/HER2-negative disease. Preserve the source's clinical
definitions and record discordant, equivocal and missing receptors. Exclude these from
the primary contrast rather than imputing labels from the transcriptome. Gene-expression
subtypes such as PAM50 are a separate variable, not substitutes for receptor status.

The pseudobulk estimand is the subtype-associated log2 expression ratio in a reviewed
cell type, adjusted for the declared covariates. The bulk estimand is the total-tissue
ratio; it can reflect cellular composition as well as expression within cells. Effect
agreement is triangulation, not evidence that the two estimands are identical.

## Statistical units and models

| Analysis | Independent unit | Implemented model / summary |
|---|---|---|
| scRNA differential expression | Patient within one cell type | Sum raw UMI counts; PyDESeq2 negative-binomial GLM |
| Bulk differential expression | One selected tumor sample per patient | PyDESeq2 negative-binomial GLM |
| Spatial program pattern | Section, exploratory within-section statistic | Gap-aware symmetric neighbors; Moran's I and label permutations |
| Spatial compartment difference | Patient after equal section averaging | Tumor minus stroma relative T-cell mixture; percentile bootstrap |
| Cross-assay effects | Same frozen gene ID within declared contrast | Effect tables and descriptive plots; no pooled-assay p values |

Default formula: `~ batch + condition`; contrast `[condition, TNBC, HRpos_HER2neg]`.
This is a minimal template, not a guarantee against confounding. Age, stage, center,
purity and technology require prespecified covariate review. Covariates used in the
formula must be complete and have appropriate types. Rank-deficient designs and
insufficient residual degrees of freedom stop the run. Between-patient subtype models
must not add patient fixed effects, which would alias the condition.

Technical replicates are summed at patient/type level only when condition and batch
are consistent. A patient spanning batches stops this engine: use a deliberately
specified repeated-measures method rather than silently assigning a batch. The bulk
engine rejects duplicate patient IDs; paired tumor-normal/longitudinal designs need
a separate patient-blocked model and are not claimed as supported here.

## Filtering and normalization

Retain raw integer counts for inference. Library-size normalization and log1p are
for embeddings and descriptive plots. Do not use integrated embeddings, scaled genes,
TPM, FPKM, ComBat outputs or log counts as negative-binomial response data.

Capture QC and pseudobulk cell thresholds are recorded in YAML. Default real-run
pseudobulks need >=30 cells and >=4 patients in each group for an eligible type.
The minimum is an execution guard, not proof of statistical power. Gene filtering
requires count >=10 in >=4 samples by default. Missing types remain missing; no
zero-count pseudobulk patient is invented. Patients without a detected type change
the target population; report this selection and assess sensitivity.

PyDESeq2 estimates size factors and dispersions, refits Cook's outliers when eligible,
and computes Wald tests. Preserve NA p values due to filtering/outliers. Effects are
**unshrunk** and labeled accordingly. For manuscript emphasis on top gene effects,
add a validated shrinkage analysis and concordance check with the target DESeq2 version.

## Multiplicity and effect reporting

The main pseudobulk exploratory family contains every tested gene × cell-type contrast.
`q_global` applies BH to all finite Wald p values in that family. `padj` is also retained
as PyDESeq2's within-type adjusted value and must not replace `q_global` in primary
multi-type statements. Bulk genes constitute a separate family. Moran tests across
section × program form another explicitly labeled exploratory family.

Report estimate, standard error, CI where available, adjusted p, independent patient n,
genes/cells included and exclusions. Figure highlighting uses q<0.05 and |log2FC|>1;
this is a display rule, not a claim that all important biology exceeds that cutoff.
The prototype marker-set ORA uses tested genes as background, signed hits, and BH
across all source × program × direction tests in `q_all_program_tests`. The four small
author-defined marker sets are not comprehensive pathways or validated biomarkers.

## Spatial methods and limitations

Build each spatial graph within one section. Connect k nearest spots up to 1.8 times
the median nearest-neighbor distance; symmetrize and remove self loops. Inspect tissue
gaps and use array topology/physical distance sensitivity when appropriate. Pixel-space
distances are not micrometers without a verified calibration. No edges cross slides.

Moran's I permutation tests use exchangeable spot labels within the section. This tests
spatial randomness, not a mechanistic ligand-receptor process. Compartment structure,
spatial autocorrelation in covariates and tissue geometry can violate a naive null.
Confirmatory colocalization would need a justified restricted/block spatial null and
patient replication; the present workflow intentionally does not claim that inference.

The NNLS CPU baseline fits shared-gene CPM signatures, averages reference patients
equally within type and excludes overlapping target patients. It reports relative
**RNA contributions**, reconstruction residuals and matrix condition numbers. It has
no calibrated posterior and is not the default for real-run mapping. Full-array dense
NNLS can be memory intensive; use the GPU workflow for larger cohorts.

The optional cell2location adapter exports mean, 5th and 95th percentile abundances,
model checkpoints and optimization histories. Only posterior means are normalized
for relative-composition displays. Quantiles are not normalized and relabeled as
fractions. Review histology-informed cells-per-spot, hyperparameter sensitivity,
residuals, convergence and reference completeness. Abundance uncertainty is not fully
propagated through the current patient bootstrap; report that limitation.

For compartment summaries, preserve ambiguous/mixed source annotations separately.
Average sections equally within each patient's compartment and compute one difference.
Bootstrap independent patients (2,000 draws in the real template); fewer than 3 patients
produce no interval. Do not interpret a relative-mixture difference as absolute
infiltration density. Orthogonal cells/mm² measurements define the proposed primary
experimental endpoint and require a separately implemented prespecified model.

## Validation, leakage and sensitivity

Keep discovery and validation patient sets disjoint. Freeze gene sets, thresholds,
hyperparameters and effect direction before external evaluation. Matching public
scRNA/spatial specimens are not independent replication. TCGA patients constitute
a separate cohort, but total-tissue expression alone cannot establish spatial biology.

Required real-study sensitivity analyses: donor leave-one-out effects; QC threshold
grid fixed before results; batch/center-aware design; subtype definitions; cell-label
review; doublet/ambient-RNA treatment; reference donor folds; spatial graph radius;
cells-per-spot priors; compartment definitions; purity effects and sampling bias.
These are requirements for a completed research study, not a claim that all were run.
The repository supplies a configuration-sweep runner for QC/graph exploratory checks.

No survival signature, optimum cutpoint, drug-response prediction, trajectory direction
or ligand-receptor causality is fabricated to increase figure count. Add such modules
only for a defined endpoint, appropriate independent validation and adequate sample size.

Sources: [DE best practices](https://www.sc-best-practices.org/conditions/differential-gene-expression/),
[PyDESeq2 workflow](https://pydeseq2.readthedocs.io/en/stable/auto_examples/plot_minimal_pydeseq2_pipeline.html),
[cell2location tutorial](https://cell2location.readthedocs.io/en/latest/notebooks/cell2location_tutorial.html).

