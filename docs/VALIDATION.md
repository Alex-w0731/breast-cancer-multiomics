# Validation record / 验证记录

Date: 2026-09-06. Overall assessment: **Share with caveats — validated software
demonstration, not a validated real-patient scientific conclusion.**

## Executed checks

| Check | Observed result |
|---|---|
| Package installation | Editable installation with dev and workflow extras completed in isolated Python 3.12.14 on Windows |
| Unit/statistical tests | 25 passed; 17 dependency deprecation/future warnings retained in the log |
| Static checks | Ruff passed for src, tests and scripts |
| CPU end-to-end workflow | demo-data → singlecell → pseudobulk → bulk → spatial → integrate → figures completed |
| Demonstration size | 12 synthetic patients; 1,920 cells; 768 spots; 12 bulk samples; 800 genes |
| Figure exports | 46 distinct figures × PNG/PDF/SVG = 138 main exports, plus web previews and contact sheet |
| Public spatial metadata download | 217,533 bytes; publisher MD5 matched; six source metadata tables inspected |
| GDC discovery | Live open TCGA-BRCA STAR-count discovery and pagination completed; no full expression cohort downloaded |
| Snakemake | DAG dry-run completed locally; Linux scheduled execution is delegated to the included GitHub CI |

The final source and installed package versions are recorded alongside the audit
JSON. `requirements-tested-windows.txt` is an environment inventory for the tested
Windows runtime, not a hash-locked cross-platform environment. Core dependencies
are explicitly pinned in pyproject.toml. Linux CI installs these pins independently.

## What the tests verify

- Noninteger, negative, nonfinite and empty count matrices are rejected.
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
full-size views are included in the gallery. All empirical-looking demo plots carry
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

