"""Discover TCGA-BRCA STAR count files, then fetch an explicitly curated selection.

Discovery never fabricates receptor status or picks an aliquot silently. Selection
TSV must contain file_id, sample_id, patient_id, condition, batch and optional
prespecified covariates. Use the source file-to-case mapping for verification.
"""

import argparse
import json

import pandas as pd

from bcflow.core import mkdir, write_json
from bcflow.download import download, session


def discover(out):
    filters = {"op": "and", "content": [{"op": "in", "content": {"field": key, "value": val}}
        for key, val in [
            ("cases.project.project_id", ["TCGA-BRCA"]),
            ("data_type", ["Gene Expression Quantification"]),
            ("analysis.workflow_type", ["STAR - Counts"]),
            ("access", ["open"]),
            ("cases.samples.sample_type", ["Primary Tumor"]),
        ]]}
    fields = "file_id,file_name,md5sum,file_size,cases.submitter_id,cases.samples.submitter_id,cases.samples.sample_type"
    offset, hits = 0, []
    s = session()
    while True:
        response = s.get("https://api.gdc.cancer.gov/files", params={
            "filters": json.dumps(filters), "fields": fields, "format": "JSON",
            "size": 500, "from": offset, "sort": "file_id:asc"}, timeout=(15, 120))
        response.raise_for_status()
        data = response.json()["data"]
        hits.extend(data["hits"])
        offset += len(data["hits"])
        if offset >= data["pagination"]["total"]:
            break
        if not data["hits"]:
            raise ValueError("Incomplete GDC pagination")
    write_json(out / "files_snapshot.json", hits)
    write_json(out / "query.json", {"filters": filters, "fields": fields})
    status = s.get("https://api.gdc.cancer.gov/status", timeout=60)
    status.raise_for_status()
    write_json(out / "gdc_release.json", status.json())
    rows = []
    for hit in hits:
        for case in hit.get("cases", []):
            rows.append({"file_id": hit["file_id"], "file_name": hit["file_name"],
                         "patient_id": case["submitter_id"], "md5": hit["md5sum"],
                         "bytes": hit["file_size"],
                         "sample_ids": ";".join(x["submitter_id"] for x in case.get("samples", [])),
                         "sample_types": ";".join(x["sample_type"] for x in case.get("samples", []))})
    pd.DataFrame(rows).to_csv(out / "candidates.tsv", sep="\t", index=False)


def fetch(out, selection, destination):
    selected = pd.read_csv(selection, sep="\t", dtype={"file_id": str})
    required = {"file_id", "sample_id", "patient_id", "condition", "batch"}
    if not required.issubset(selected) or selected[list(required)].isna().any().any():
        raise ValueError("Curated selection is missing required fields")
    if selected.patient_id.duplicated().any() or selected.file_id.duplicated().any():
        raise ValueError("Select exactly one file per patient for the between-patient model")
    if selected.sample_id.duplicated().any():
        raise ValueError("Duplicate sample IDs")
    candidates = pd.read_csv(out / "candidates.tsv", sep="\t").set_index("file_id")
    matrices, genes = [], None
    for row in selected.itertuples(index=False):
        item = candidates.loc[row.file_id]
        if isinstance(item, pd.DataFrame) or item.patient_id != row.patient_id:
            raise ValueError("Ambiguous file-to-patient mapping")
        if row.sample_id not in str(item.sample_ids).split(";"):
            raise ValueError("Selected sample_id does not match the GDC file metadata")
        path = download(f"https://api.gdc.cancer.gov/data/{row.file_id}",
                        out / "counts" / f"{row.file_id}.tsv", item.md5,
                        max_bytes=int(item.bytes) + 1024)
        table = pd.read_csv(path, sep="\t", comment="#")
        if not {"gene_id", "gene_name", "unstranded"}.issubset(table):
            raise ValueError("Unexpected GDC STAR-count schema; use raw unstranded column")
        table = table[table.gene_id.str.startswith("ENSG")].copy()
        table["gene_id"] = table.gene_id.str.replace(r"\.\d+$", "", regex=True)
        if table.gene_id.duplicated().any():
            raise ValueError("Gene ID normalization introduced duplicate features")
        table = table.set_index("gene_id")
        if genes is None:
            genes = table[["gene_name"]].rename(columns={"gene_name": "gene_symbol"})
        elif not genes.index.equals(table.index) or not genes.gene_symbol.equals(table.gene_name):
            raise ValueError("Inconsistent gene annotation across GDC files")
        matrices.append(table.unstranded.rename(row.sample_id))
    if not matrices:
        raise ValueError("Selection is empty")
    dest = mkdir(destination)
    pd.DataFrame(matrices).to_csv(dest / "bulk_counts.csv")
    selected.set_index("sample_id").to_csv(dest / "bulk_metadata.csv")
    genes.to_csv(dest / "bulk_genes.csv")


def main():
    p = argparse.ArgumentParser()
    p.add_argument("command", choices=["discover", "fetch"])
    p.add_argument("--output", default="data/gdc")
    p.add_argument("--selection")
    p.add_argument("--destination", default="data/processed")
    args = p.parse_args()
    out = mkdir(args.output)
    if args.command == "discover":
        discover(out)
    elif args.selection:
        fetch(out, args.selection, args.destination)
    else:
        p.error("fetch requires --selection")


if __name__ == "__main__":
    main()
