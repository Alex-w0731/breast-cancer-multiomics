# Validation record / 验证记录

Date: 2026-09-06. Overall assessment: **Share with caveats — validated software
demonstration, not a validated real-patient scientific conclusion.**

## Executed checks

| Check | Observed result |
|---|---|
| Package installation | Editable installation with dev and workflow extras completed in isolated Python 3.12.14 on Windows |
| Unit/statistical tests | 28 passed; 17 dependency deprecation/future warnings retained in the log |
| Static checks | Ruff passed for src, tests and scripts |
| CPU end-to-end workflow | demo-data → singlecell → pseudobulk → bulk → spatial → integrate → figures completed |
| Demonstration size | 12 synthetic patients; 1,920 cells; 768 spots; 12 bulk samples; 800 genes |
| Figure exports | 46 distinct figures × PNG/PDF/SVG = 138 main exports, plus web previews and contact sheet |
| Public spatial metadata download | 217,533 bytes; publisher MD5 matched; six source metadata tables inspected |
| GDC discovery | Live open TCGA-BRCA STAR-count discovery returned 1,111 candidate files from 1,095 distinct patients; no full expression cohort downloaded |
| Snakemake | DAG dry-run completed locally; the complete CPU DAG also executed successfully in GitHub CI on Ubuntu, with two cores |

The final code passed installation, lint, all 28 tests, the complete Snakemake DAG
and artifact upload: [GitHub Actions run 34005821215](https://github.com/Alex-w0731/breast-cancer-multiomics/actions/runs/34005821215),
commit `03d2e1a4f1572f87c128b9f86d53f12973aed06a`. The full observed job/step
summary is archived in [validation/github_ci.json](validation/github_ci.json).
Download the full-resolution figures and audit files from that run's artifact;
GitHub retains this configured CI artifact for 14 days. The supplied local ZIP
preserves the exports, synthetic inputs, source tables and a SHA-256 file manifest.
GDC discovery recorded Data Release 46.0 (August 10, 2026); candidate counts are a
dated API snapshot, not the number of clinically eligible independent samples.

The final source and installed package versions are recorded alongside the audit
JSON. `requirements-tested-windows.txt` is an environment inventory for the tested
Windows runtime, not a hash-locked cross-platform environment. Core dependencies
are explicitly pinned in pyproject.toml. Linux CI installs these pins independently.

## What the tests verify

- Noninteger, negative, nonfinite and empty count matrices are rejected.
- Duplicate CSV gene headers and missing patient IDs are rejected before model fitting;
  nonfinite audit values serialize as JSON null instead of invalid NaN literals.
- Synthetic AnnData is blocked from a real run.
- Pseudobulk uses raw sums and conserves counts; technical/cell replication is not
  mistaken for patient replication; ambiguous cross-batch patients are rejected.
- Confounded/rank-deficient designs, duplicate patients and misordered metadata stop.
- A real negative-binomial test recovers the planted direction and significant
  effects in an independent simulated fixture.
- Global BH adjustment preserves missing tests; gene-set background excludes untested genes.
- Spatial graphs have no self loops or long bridges across separated tissue islands.
- Moran's I matches the directly computed dense formula; constant fields return no fake test.
- NNLS recovers a known pure mixture and excludes the target donor from its reference.
- Bootstrap does not invent uncertainty for two donors; unsafe archive paths/links are rejected.

## Visual review

The exported contact sheet was inspected for missing figures, illegible labels,
incorrect coordinate orientation, inconsistent status labels and clipping. Individual
full-size views were reviewed locally; the repository gallery uses browsing previews.
The marker dot plot includes a detection-fraction size legend and an observed color
range; patient intervals are drawn from their actual endpoints. All empirical-looking demo plots carry
SYNTHETIC DATA - SOFTWARE DEMO ONLY. The power plot is explicitly hypothetical.
No image-generation model was used to invent histology or biological measurements.

## Issues fixed during validation

The Scanpy igraph Leiden backend produced repeated Windows integer-range errors.
The code now explicitly uses the tested leidenalg backend and fixed UMAP/Numba/
PyNNDescent versions. The corrected end-to-end log contains no integer-range error.
Source/config dependencies and all consumed core tables are declared in Snakemake;
old patient-compartment outputs are removed before writing the current integration.
Application caches are directed to task-local work directories where possible.

Remaining warnings come from dependencies (including PyDESeq2 dtype deprecations
and occasional dispersion-trend fallback); they are preserved, not silently hidden.
Real-data fits must inspect convergence, model diagnostics and sensitivity rather
than treating successful process exit as sufficient evidence.

## Unrun work and research limitations

- No complete GSE176078/TCGA/Visium expression cohort analysis was completed.
- Full raw FASTQ alignment, sample clinical curation, receptor reconciliation and
  GDC-selected expression downloads remain real-study steps.
- GPU cell2location fitting was not executed; its adapter is source-reviewed but
  requires a separately tested environment and posterior diagnostics.
- Production Scrublet was not benchmarked on this public cohort.
- No independent experimental cohort, histology endpoint model, survival model,
  causal perturbation or prospective study outcome is claimed.
- NNLS is an exploratory RNA-mixture baseline. Spatial bootstrap does not propagate
  the complete deconvolution posterior uncertainty.
- Public spatial metadata contain only four TNBC and two ER-labeled samples;
  ER labels and mixed pathology require further review, not automatic relabeling.

These items block a publication claim about real biology, not distribution of the
clearly labeled software and proposed experimental protocol. See the complete
[publication checklist](publication_checklist.md) and [SAP](statistical_analysis_plan.md).
