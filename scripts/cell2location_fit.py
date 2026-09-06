"""Optional GPU mapping adapter. Requires a separately validated cell2location environment.

Use disjoint reference patients for validation. This expensive optional adapter
is source-reviewed, but is not included in the CPU demonstration's execution claim.
"""

import argparse
import json
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--reference", required=True)
    p.add_argument("--spatial", required=True)
    p.add_argument("--output-dir", required=True)
    p.add_argument("--cells-per-spot", type=float, required=True,
                   help="Estimate from nuclei counts in the corresponding H&E section")
    p.add_argument("--reference-epochs", type=int, default=250)
    p.add_argument("--spatial-epochs", type=int, default=30000)
    p.add_argument("--accelerator", choices=["cpu", "gpu"], default="gpu")
    p.add_argument("--seed", type=int, default=7301)
    args = p.parse_args()
    import cell2location
    import scvi
    from cell2location.models import Cell2location, RegressionModel
    from cell2location.utils.filtering import filter_genes
    from importlib.metadata import version

    scvi.settings.seed = args.seed
    ref, st = ad.read_h5ad(args.reference), ad.read_h5ad(args.spatial)
    overlap = set(ref.obs.patient_id) & set(st.obs.patient_id)
    if overlap:
        raise ValueError("Reference and target patients overlap. Use donor-disjoint folds for validation")
    if args.cells_per_spot <= 0:
        raise ValueError("cells-per-spot must be a positive histology-informed estimate")
    for a in [ref, st]:
        if "counts" not in a.layers:
            raise ValueError("Raw counts layer required")
        a.X = a.layers["counts"].copy()
        v = a.X.data if hasattr(a.X, "toarray") else np.asarray(a.X)
        if not np.isfinite(v).all() or (v < 0).any() or not (v == np.floor(v)).all():
            raise ValueError("Raw integer counts required")
    out = Path(args.output_dir)
    out.mkdir(parents=True, exist_ok=True)
    selected = filter_genes(ref, cell_count_cutoff=5, cell_percentage_cutoff2=.03,
                            nonz_mean_cutoff=1.12)
    genes = ref[:, selected].var_names.intersection(st.var_names, sort=False)
    if len(genes) < 500:
        raise ValueError("Fewer than 500 shared genes")
    ref, st = ref[:, genes].copy(), st[:, genes].copy()
    RegressionModel.setup_anndata(ref, batch_key="sample_id", labels_key="cell_type")
    model_ref = RegressionModel(ref)
    model_ref.train(max_epochs=args.reference_epochs, accelerator=args.accelerator,
                     batch_size=2500, train_size=1)
    ref = model_ref.export_posterior(ref, sample_kwargs={"num_samples": 1000,
                                     "batch_size": 2500, "accelerator": args.accelerator})
    names = list(ref.uns["mod"]["factor_names"])
    columns = [f"means_per_cluster_mu_fg_{x}" for x in names]
    slot = ref.varm["means_per_cluster_mu_fg"] if "means_per_cluster_mu_fg" in ref.varm else ref.var
    signatures = pd.DataFrame(slot[columns], index=ref.var_names)
    signatures.columns = names
    signatures.to_csv(out / "reference_signatures.csv")
    model_ref.save(str(out / "reference_model"), overwrite=False)
    ref.write_h5ad(out / "reference_posterior.h5ad")
    Cell2location.setup_anndata(st, batch_key="section_id")
    model = Cell2location(st, cell_state_df=signatures,
                          N_cells_per_location=args.cells_per_spot, detection_alpha=20)
    model.train(max_epochs=args.spatial_epochs, batch_size=None, train_size=1,
                 accelerator=args.accelerator)
    st = model.export_posterior(st, sample_kwargs={"num_samples": 1000,
                               "batch_size": st.n_obs, "accelerator": args.accelerator})
    model.save(str(out / "spatial_model"), overwrite=False)
    st.write_h5ad(out / "spatial_posterior.h5ad")
    for statistic in ["means", "q05", "q95"]:
        values = st.obsm[f"{statistic}_cell_abundance_w_sf"]
        frame = pd.DataFrame(np.asarray(values), index=st.obs_names,
                             columns=st.uns["mod"]["factor_names"])
        frame.to_csv(out / f"abundance_{statistic}.csv")
    for name, model_object in [("reference", model_ref), ("spatial", model)]:
        for history, value in model_object.history.items():
            pd.DataFrame(value).to_csv(out / f"{name}_{history}.csv")
    (out / "run_parameters.json").write_text(json.dumps({**vars(args),
        "cell2location": version("cell2location"), "scvi-tools": version("scvi-tools"),
        "shared_genes": len(genes), "package_path": cell2location.__file__}, indent=2))


if __name__ == "__main__":
    main()
