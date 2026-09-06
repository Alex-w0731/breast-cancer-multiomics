"""Convert explicit matrix/feature/metadata contracts into a validated AnnData.

Matrix Market input must be genes x observations. features.tsv has header
gene_id,gene_symbol; metadata.csv uses observation ID as its first column.
Spatial metadata also requires pixel_x,pixel_y from tissue_positions.csv.
This adapter intentionally does not guess patient IDs or clinical subtype.
"""

import argparse

import anndata as ad
import pandas as pd
from scipy import io, sparse

from bcflow.core import mkdir, sha256, validate_adata


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--matrix", required=True)
    p.add_argument("--features", required=True)
    p.add_argument("--barcodes", required=True)
    p.add_argument("--metadata", required=True)
    p.add_argument("--modality", choices=["singlecell", "spatial"], required=True)
    p.add_argument("--output", required=True)
    args = p.parse_args()
    genes = pd.read_csv(args.features, sep="\t", index_col="gene_id")
    barcodes = pd.read_csv(args.barcodes, sep="\t", header=None, dtype=str)[0]
    meta = pd.read_csv(args.metadata, index_col=0)
    if not meta.index.is_unique or not barcodes.is_unique or not genes.index.is_unique:
        raise ValueError("Duplicate IDs in inputs")
    if set(meta.index) != set(barcodes):
        raise ValueError("Matrix barcodes and metadata must match exactly")
    x = sparse.csr_matrix(io.mmread(args.matrix)).T.tocsr()
    if x.shape != (len(barcodes), len(genes)):
        raise ValueError("Matrix orientation/shape differs from barcode and feature files")
    a = ad.AnnData(x, obs=meta.loc[barcodes].copy(), var=genes)
    a.layers["counts"] = a.X.copy()
    a.uns["data_status"] = "PUBLIC_OR_CONSENTED_REAL"
    a.uns["input_sha256"] = {key: sha256(getattr(args, key)) for key in
                            ["matrix", "features", "barcodes", "metadata"]}
    if args.modality == "spatial":
        a.obsm["spatial"] = a.obs[["pixel_x", "pixel_y"]].to_numpy(dtype=float)
    validate_adata(a, args.modality, "real")
    from pathlib import Path
    mkdir(Path(args.output).parent)
    a.write_h5ad(args.output, compression="gzip")


if __name__ == "__main__":
    main()
