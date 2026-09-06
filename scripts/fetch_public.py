"""Download selected Wu spatial files with publisher MD5; GEO accession metadata."""

import argparse
from pathlib import Path

from bcflow.download import download, safe_extract

SPATIAL = {
    "metadata.tar.gz": "1d0b34ead70635c094ac74ac58a88d68",
    "filtered_count_matrices.tar.gz": "ea1220f0606c4d4e9307468fb2d5427a",
    "spatial.tar.gz": "29b02ae433d459d19dc52b75a529e8e9",
}


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--part", choices=["metadata", "spatial", "singlecell"], required=True)
    p.add_argument("--output", default="data/public")
    p.add_argument("--extract", action="store_true")
    args = p.parse_args()
    out = Path(args.output)
    if args.part == "singlecell":
        # GEO does not advertise a publisher MD5 here; save a local SHA-256 receipt.
        filename = "GSE176078_Wu_etal_2021_BRCA_scRNASeq.tar.gz"
        path = download("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE176nnn/GSE176078/suppl/"
                        + filename, out / filename, max_bytes=1_000_000_000)
        if args.extract:
            safe_extract(path, out / "singlecell")
    else:
        selected = ["metadata.tar.gz"] if args.part == "metadata" else list(SPATIAL)
        for filename in selected:
            path = download(f"https://zenodo.org/records/4739739/files/{filename}?download=1",
                            out / filename, SPATIAL[filename], max_bytes=500_000_000)
            if args.extract:
                safe_extract(path, out / filename.removesuffix(".tar.gz"))


if __name__ == "__main__":
    main()

