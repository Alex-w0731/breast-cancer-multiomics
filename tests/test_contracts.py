import numpy as np
import pandas as pd
import pytest
from scipy import sparse
import anndata as ad

from bcflow.core import bh, validate_counts, validate_adata
from bcflow.differential import aggregate, check_design
from bcflow.integration import bootstrap_mean, ora
from bcflow.spatial import graph, moran


@pytest.mark.parametrize("bad", [np.array([[1.1, 2]]), np.array([[-1, 2]]),
                                  np.array([[np.nan, 2]]), np.array([[0, 0]])])
def test_reject_invalid_counts(bad):
    with pytest.raises(ValueError):
        validate_counts(bad)


def test_raw_counts_accepted_dense_and_sparse():
    for matrix in [np.array([[0, 2], [3, 4]]), sparse.csr_matrix([[0, 2], [3, 4]])]:
        validate_counts(matrix)


def test_global_bh_handles_missing():
    assert np.allclose(bh([.01, .04, .03, np.nan])[:3], [.03, .04, .04])
    assert np.isnan(bh([np.nan])[0])


def example():
    obs = pd.DataFrame({"patient_id": ["p1", "p1", "p2"], "cell_type": ["T", "T", "T"],
                        "sample_id": ["s1", "s1", "s2"], "condition": ["A", "A", "B"],
                        "batch": ["b1", "b1", "b1"], "annotation_source": "reviewed"},
                       index=["c1", "c2", "c3"])
    a = ad.AnnData(sparse.csr_matrix([[1, 2], [3, 4], [5, 6]]), obs=obs,
                   var=pd.DataFrame({"gene_symbol": ["A", "B"]}, index=["A", "B"]))
    a.layers["counts"] = a.X.copy()
    return a


def test_pseudobulk_sums_not_averages():
    counts, meta, dropped = aggregate(example(), 1)
    assert counts.loc["p1::T"].tolist() == [4, 6]
    assert counts.to_numpy().sum() == 21
    assert meta.loc["p1::T", "n_cells"] == 2
    assert len(meta) == 2
    assert dropped.empty


def test_patient_cross_batch_is_not_silently_collapsed():
    a = example()
    a.obs.loc["c2", "batch"] = "b2"
    with pytest.raises(ValueError, match="multiple batch"):
        aggregate(a, 1)


def test_synthetic_cannot_enter_real():
    a = example()
    a.uns["data_status"] = "SYNTHETIC"
    with pytest.raises(ValueError, match="Synthetic"):
        validate_adata(a, "singlecell", "real")


def design_fixture():
    md = pd.DataFrame({"patient_id": [f"p{i}" for i in range(8)],
                       "condition": ["A", "B"] * 4, "batch": ["x", "x", "y", "y"] * 2},
                      index=[f"s{i}" for i in range(8)])
    counts = pd.DataFrame(np.ones((8, 30), dtype=int), index=md.index)
    cfg = {"contrast": ["condition", "B", "A"], "design": "~ batch + condition",
           "analysis": {"min_donors_per_group": 3}}
    return counts, md, cfg


def test_balanced_design_valid():
    counts, md, cfg = design_fixture()
    check_design(counts, md, cfg)


def test_confounded_design_rejected():
    counts, md, cfg = design_fixture()
    md["batch"] = md.condition
    with pytest.raises(ValueError, match="Rank-deficient"):
        check_design(counts, md, cfg)


def test_pseudoreplication_rejected():
    counts, md, cfg = design_fixture()
    md.loc["s1", "patient_id"] = "p0"
    with pytest.raises(ValueError, match="independent observation"):
        check_design(counts, md, cfg)


def test_order_mismatch_rejected():
    counts, md, cfg = design_fixture()
    with pytest.raises(ValueError, match="ordered sample"):
        check_design(counts, md.iloc[::-1], cfg)


def test_graph_does_not_bridge_distant_tissue():
    coords = np.array([[0, 0], [1, 0], [0, 1], [1, 1],
                        [100, 100], [101, 100], [100, 101], [101, 101]])
    w = graph(coords)
    assert w[:4, 4:].nnz == 0
    assert (w != w.T).nnz == 0
    assert w.diagonal().sum() == 0


def test_moran_agrees_with_dense_formula():
    coords = np.array([[i, j] for i in range(5) for j in range(5)])
    w = graph(coords)
    values = coords[:, 0].astype(float)
    value, p = moran(values, w, 99, np.random.default_rng(3))
    z = values - values.mean()
    expected = len(z) / w.sum() * (z @ w.toarray() @ z) / (z @ z)
    assert np.isclose(value, expected)
    assert .01 <= p <= 1
    assert value > 0


def test_constant_spatial_program_has_no_test():
    w = graph(np.array([[0, 0], [1, 0], [0, 1], [1, 1]]))
    assert np.isnan(moran(np.ones(4), w, 99, np.random.default_rng(1))[0])


def test_bootstrap_no_fake_interval_for_two_donors():
    _, low, high = bootstrap_mean([1, 2], np.random.default_rng(1))
    assert np.isnan(low) and np.isnan(high)


def test_ora_uses_tested_background():
    df = pd.DataFrame({"gene_symbol": ["COL1A1", "COL1A2", "X", "Y", "FN1"],
                       "pvalue": [.001, .001, .6, .7, np.nan], "padj": [.01, .01, .9, .9, np.nan],
                       "log2FoldChange": [2, 2, 0, 0, 2]})
    result = ora(df)
    row = result.query("program == 'ECM' and direction == 'up'").iloc[0]
    assert row.universe_size == 4
    assert row.set_size_tested == 2
    assert row.overlap == 2


def test_strict_json_preserves_missing_and_numpy_types(tmp_path):
    import json
    from bcflow.core import write_json
    path = tmp_path / "audit.json"
    write_json(path, {"ci": np.nan, "n": np.int64(3), "valid": np.bool_(True)})
    content = path.read_text()
    assert "NaN" not in content
    assert json.loads(content) == {"ci": None, "n": 3, "valid": True}


def test_duplicate_gene_csv_header_rejected(tmp_path):
    from bcflow.core import read_count_csv
    path = tmp_path / "counts.csv"
    path.write_text(",geneA,geneA\nsample1,10,20\n")
    with pytest.raises(ValueError, match="Duplicate"):
        read_count_csv(path)


def test_missing_patient_id_rejected():
    counts, md, cfg = design_fixture()
    md.loc["s0", "patient_id"] = None
    with pytest.raises(ValueError, match="patient_id"):
        check_design(counts, md, cfg)
