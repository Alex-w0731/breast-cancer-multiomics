# Data sources and evidence register

Accessed/checked: **2026-09-06**. Availability may change; archive each actual download
and metadata snapshot with a checksum and license/access notes.

| Layer | Source | Intended use | Important limitation |
|---|---|---|---|
| scRNA-seq | [GSE176078](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE176078) | Reviewed breast tumor cell reference and patient pseudobulks | Public processed data; original raw human reads have controlled access |
| Spatial | [Zenodo 4739739](https://zenodo.org/records/4739739), DOI 10.5281/zenodo.4739739 | Six Visium tumors, positions and pathology metadata | Same study; overlapping patients are not independent replication |
| bulk | [GDC TCGA-BRCA](https://portal.gdc.cancer.gov/projects/TCGA-BRCA) | Independent bulk expression cohort | Clinical receptor mapping and aliquot choice require review |
| Original study | [Wu et al. 2021](https://doi.org/10.1038/s41588-021-00911-1) | Biological/data context | A published dataset does not validate a new hypothesis automatically |
| Raw data access | [EGA EGAS00001005173](https://ega-archive.org/studies/EGAS00001005173) | Optional authorized raw-read reprocessing | Controlled access; no access restrictions bypassed |

## Verified real-data smoke check

Downloaded only the 217,533-byte spatial metadata archive for the present software
validation. Publisher MD5 `1d0b34ead70635c094ac74ac58a88d68` matched.
Local SHA-256: `beb8c05063f279fef805349ee8cc0714957576fbfb4af7c119dc46ec5cc8bea3`.
The count matrices and whole TCGA cohort were not analyzed in this delivery.

| Metadata file | Rows in archived table | Stored subtype |
|---|---:|---|
| 1142243F_metadata.csv | 4,784 | TNBC |
| 1160920F_metadata.csv | 4,895 | TNBC |
| CID4290_metadata.csv | 2,432 | ER |
| CID4465_metadata.csv | 1,211 | TNBC |
| CID44971_metadata.csv | 1,162 | TNBC |
| CID4535_metadata.csv | 1,127 | ER |

Total: 15,611 metadata rows in six files. These are archive table counts, not the number
of spots that will pass this pipeline's QC. The `ER` label alone does not establish
HER2 negativity. Several pathology labels are mixed tumor/stroma/lymphocyte regions;
mapping them to pure compartments without review would change the estimand.

## Data and tool references

- [GDC query API](https://docs.gdc.cancer.gov/API/Users_Guide/Search_and_Retrieval/):
  file discovery, identifiers and metadata snapshots.
- [GDC STAR Counts](https://docs.gdc.cancer.gov/Encyclopedia/pages/STAR_Counts/):
  harmonized RNA expression count output.
- [PyDESeq2 API](https://pydeseq2.readthedocs.io/en/stable/api/docstrings/pydeseq2.dds.DeseqDataSet.html):
  counts, metadata and design formula interface.
- [Single-cell differential-expression best practices](https://www.sc-best-practices.org/conditions/differential-gene-expression/):
  patient-level replication and pseudobulk reasoning.
- [cell2location tutorial](https://cell2location.readthedocs.io/en/latest/notebooks/cell2location_tutorial.html):
  reference NB regression and spatial posterior modeling.
- [Scanpy Leiden API](https://scanpy.readthedocs.io/en/stable/api/scanpy.tl.leiden.html):
  graph clustering. This repository explicitly uses the leidenalg backend.
- [Scanpy Windows backend issue](https://github.com/scverse/scanpy/issues/3028):
  reason the igraph backend is not used in the tested environment.

Code is newly written for this repository. Public records retain their original
ownership and restrictions; the repository's MIT software license does not relicense data.
Do not upload identifiable clinical metadata, FASTQ/BAM or access tokens to GitHub.

