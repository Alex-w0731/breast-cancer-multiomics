# Raw sequencing preprocessing boundary

The tested BCFlow entry point is raw **gene-count matrices**, not FASTQ.
Below are reviewable upstream commands for compatible assays. No FASTQ run,
alignment benchmark or commercial 10x installation was completed in this delivery.

For a real study, freeze GRCh38 FASTA and GTF releases, transcript filtering,
gene-ID mapping and SHA-256 manifests. The scRNA, spatial and bulk branches need
compatible gene annotation; a matching genome name alone is insufficient.

## Single-cell 10x gene expression

For a compatible singleplex 3'/5' GEX library:

```bash
cellranger count --id=SC_SAMPLE_ID \
  --transcriptome=/reference/FROZEN_GRCH38_10X \
  --fastqs=/data/fastqs/SC_SAMPLE_ID \
  --sample=FASTQ_SAMPLE_PREFIX --create-bam=true \
  --localcores=16 --localmem=64
```

Record Cell Ranger version, chemistry, include-introns behavior, sequencing depth,
valid barcodes, mapping rate, estimated cells, saturation and ambient RNA assessment.
The command is not appropriate for every multiplexed/Flex/Multiome assay: select
the pipeline matching the actual chemistry in the [official Cell Ranger guide](https://www.10xgenomics.com/support/software/cell-ranger/latest/analysis/running-pipelines/cr-gex-count).
Retain raw and filtered feature-barcode matrices plus web_summary.html.

## Visium spatial GEX

Use the [official Space Ranger count guide](https://www.10xgenomics.com/support/software/space-ranger/latest/analysis/running-pipelines/space-ranger-count)
for the exact slide/chemistry and image type. A conventional fresh-frozen Visium
example starts from the following fields:

```bash
spaceranger count --id=SPATIAL_SAMPLE_ID \
  --transcriptome=/reference/FROZEN_GRCH38_10X \
  --fastqs=/data/fastqs/SPATIAL_SAMPLE_ID --sample=FASTQ_SAMPLE_PREFIX \
  --image=/data/images/section.tif --slide=ACTUAL_SLIDE_SERIAL \
  --area=ACTUAL_CAPTURE_AREA --create-bam=true \
  --localcores=16 --localmem=64
```

Check image alignment, tissue segmentation, fiducials and stain orientation manually.
FFPE probe-based and Visium HD assays require different preparation/parameters and
may not be comparable to the conventional Visium dataset without a defined adapter.
Keep scalefactors, tissue_positions, full-resolution image and annotated masks.

## Conventional bulk RNA-seq

Use FastQC/MultiQC for raw and processed read checks; fastp for a documented adapter/
quality-filtering rule; STAR for genome alignment with a frozen index. Example:

```bash
fastp --in1 SAMPLE_R1.fastq.gz --in2 SAMPLE_R2.fastq.gz \
  --out1 SAMPLE.clean_R1.fastq.gz --out2 SAMPLE.clean_R2.fastq.gz \
  --detect_adapter_for_pe --thread 8 --json SAMPLE.fastp.json --html SAMPLE.fastp.html
STAR --runThreadN 16 --genomeDir /reference/FROZEN_STAR_INDEX \
  --readFilesIn SAMPLE.clean_R1.fastq.gz SAMPLE.clean_R2.fastq.gz \
  --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate \
  --quantMode GeneCounts --outFileNamePrefix SAMPLE.
```

Use the correct ReadsPerGene count column for the validated library strandedness;
this must not be inferred from a filename. Remove the four special summary rows.
An alternative featureCounts approach requires correct paired-fragment and stranded
settings, with the same GTF used for interpretation. The GDC import path already uses
harmonized STAR **unstranded** values as documented by GDC; do not mix a differently
defined count matrix into that cohort without review.

Additional raw QC: contamination, uniquely/multi-mapped fraction, rRNA, gene-body
coverage, junctions, duplicate fraction, strandedness and sample swaps. FASTQ files
and clinical linkage tables stay in approved storage; publish code and permitted
deidentified derived data through the appropriate repositories.

