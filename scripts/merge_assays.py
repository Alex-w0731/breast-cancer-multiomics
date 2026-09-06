"""Merge prepared assays without silently dropping or zero-filling different gene universes."""

import argparse
from pathlib import Path

import anndata as ad

from bcflow.core import mkdir, sha256, validate_adata


def main():
    p = argparse.ArgumentParser()
    p.add_argument("inputs", nargs="+")
    p.add_argument("--modality", choices=["singlecell", "spatial"], required=True)
    p.add_argument("--output", required=True)
    args = p.parse_args()
    objects = [ad.read_h5ad(path) for path in args.inputs]
    for obj in objects:
        validate_adata(obj, args.modality, "real")
        if not objects[0].var_names.equals(obj.var_names):
            raise ValueError("Different gene universe/order: harmonize explicitly before concatenation")
        if not objects[0].var.gene_symbol.equals(obj.var.gene_symbol):
            raise ValueError("Different gene annotation")
    merged = ad.concat(objects, join="inner", merge="same", uns_merge="unique")
    validate_adata(merged, args.modality, "real")
    merged.uns["source_files_sha256"] = {path: sha256(path) for path in args.inputs}
    mkdir(Path(args.output).parent)
    merged.write_h5ad(args.output, compression="gzip")


if __name__ == "__main__":
    main()

