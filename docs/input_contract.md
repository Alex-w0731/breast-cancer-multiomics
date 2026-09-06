# Input contract

Every path in a YAML config is relative to the **repository working directory**.
The implemented comparison is between patients, with exactly two declared condition
levels. A source metadata table must be reviewed before conforming it to this contract.

## AnnData: singlecell.h5ad

- Observations are cells, variables are genes; sparse CSR recommended.
- Unique `obs_names` globally; prefix reused 10x barcodes with sample ID during preparation.
- `layers['counts']`: finite, nonnegative integer raw gene UMI counts.
- `var_names`: stable IDs from one frozen reference (prefer Ensembl gene IDs without version).
- `var['gene_symbol']`: annotation-backed gene symbols. Never derive IDs from capitalization.
- `obs`: `sample_id`, `patient_id`, `condition`, `batch`, `cell_type`, `annotation_source`.
- `annotation_source`: traceable reviewed reference / curator / revision; `automatic`,
  `unreviewed`, `unknown` and `synthetic_truth` are rejected by the real template.
- Additional donor-constant covariates may be used in the declared formula.
- Sample must map to one patient, condition and batch. Patient must have one condition.
- `uns['data_status']`: `SYNTHETIC` for demos; real data must not use this marker.

The core workflow uses reviewed labels rather than silently claiming cell annotation
is solved. First inspect the marker dot plot, source annotations and pathology. For a
new dataset without labels, establish a separate annotation pass and record its source.
Suggested major labels include Epithelial, T_cell, B_cell, Myeloid, Fibroblast,
Endothelial and other biologically justified groups. All real reference types should
be represented; the five-type synthetic fixture is not a complete breast atlas.

## AnnData: spatial.h5ad

Same counts and gene contracts; observations are spatial spots. Metadata requires
`sample_id`, `patient_id`, `condition`, `batch`, `section_id`, `in_tissue` (numeric 0/1).
`obsm['spatial']` is N×2 finite **x, y in full-resolution image pixels**.
Do not reverse the row/column coordinate order from Space Ranger.

The real template requires `pathology`, using reviewed `tumor`, `stroma` or explicit
other labels such as `mixed`, `uncertain`, `necrosis`. Mixed spots are not automatically
assigned to pure tumor. `section_id` must be globally unique. Every section maps to
one patient. Preserve `pathology_original` and the mapping version separately.

Optional histology belongs in standard AnnData `uns['spatial']` along with scalefactors;
current plotting outputs coordinate-only maps and labels them as such. They are not
claimed to be H&E overlays. Do not generate synthetic H&E to decorate real findings.

External abundance CSV: index=exact retained spot IDs, columns=reviewed cell types,
finite nonnegative posterior means, positive row sums. `T_cell` is needed for the
implemented compartment endpoint. Record method and posterior output provenance.

## Bulk tables

`bulk_counts.csv`: rows=samples, columns=stable gene IDs, first column=index.
`bulk_metadata.csv`: rows indexed by sample_id; `patient_id`, `condition`, `batch`
and all declared covariates. One sample per patient. Counts and metadata must have
the exact same sample set (the loader aligns order explicitly).
`bulk_genes.csv`: first column `gene_id`, with `gene_symbol` column for all count genes.

GDC uses the **unstranded** raw STAR-count column for the harmonized files selected here.
The script excludes special N_ summary rows, strips Ensembl version suffixes, and
rejects duplicates. It never sends TPM/FPKM into DESeq2.

## Matrix Market adapter

```bash
python scripts/prepare_matrix.py --modality singlecell \
  --matrix data/curated/counts.mtx.gz \
  --features data/curated/features.tsv \
  --barcodes data/curated/barcodes.tsv \
  --metadata data/curated/cell_metadata.csv \
  --output data/processed/singlecell.h5ad
```

Matrix is **genes × observations**. `features.tsv` has a tab-separated header with
`gene_id` and `gene_symbol`; barcode file is one column without header. Metadata CSV
first column is the matching observation ID. For spatial data add `pixel_x,pixel_y`
columns and use `--modality spatial`. Do not pass original unreviewed author metadata
without mapping its column names and clinically verifying the subtype.

For genes provided only as symbols, obtain a frozen annotation mapping, keep unambiguous
matches, and document exclusions. Never use `var_names_make_unique()` suffixes as if
they were Ensembl IDs. Differences between reference releases need an explicit
mapping table. Do not silently map ambiguous gene symbols many-to-many.

