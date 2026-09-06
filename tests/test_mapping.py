import anndata as ad
import numpy as np
import pandas as pd
import pytest
from scipy import sparse

from bcflow.spatial import deconvolve_nnls


def test_nnls_known_pure_reference_and_donor_exclusion():
    obs = pd.DataFrame({"patient_id": ["p1", "p1", "p2", "p2", "p3", "p3"],
                        "cell_type": ["A", "B"] * 3, "condition": "same", "batch": "one"},
                       index=[f"c{i}" for i in range(6)])
    x = sparse.csr_matrix([[100, 1], [1, 100]] * 3)
    ref = ad.AnnData(x, obs=obs)
    ref.layers["counts"] = x.copy()
    target = ad.AnnData(sparse.csr_matrix([[100, 1]]),
                        obs=pd.DataFrame({"patient_id": ["p1"]}, index=["s0"]))
    target.layers["counts"] = target.X.copy()
    weights, residual, details = deconvolve_nnls(ref, target, 2)
    assert weights.loc["s0", "A"] > .999
    assert np.isclose(weights.loc["s0"].sum(), 1)
    assert residual[0] < 1e-8
    assert details["excluded_overlapping_patients"] == ["p1"]


def test_duplicate_spatial_positions_rejected():
    from bcflow.spatial import graph
    with pytest.raises(ValueError, match="distinct"):
        graph([[0, 0], [0, 0], [1, 1], [2, 2]])

