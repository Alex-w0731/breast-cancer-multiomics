"""A small NB fit verifies real inference, not only mocked control flow."""

import numpy as np
import pandas as pd

from bcflow.differential import fit


def test_nb_fit_recovers_planted_direction():
    rng = np.random.default_rng(712)
    mu = np.full((16, 100), 100.0)
    mu[8:, :5] *= 5
    counts = pd.DataFrame(rng.negative_binomial(40, 40 / (40 + mu)),
                           index=[f"s{i}" for i in range(16)], columns=[f"g{i}" for i in range(100)])
    md = pd.DataFrame({"patient_id": [f"p{i}" for i in range(16)],
                       "condition": ["A"] * 8 + ["B"] * 8}, index=counts.index)
    cfg = {"design": "~ condition", "contrast": ["condition", "B", "A"],
           "analysis": {"min_donors_per_group": 3, "min_count": 10, "min_samples_gene": 3,
                        "threads": 1}}
    result, normalized = fit(counts, md, cfg)
    assert (result.loc[[f"g{i}" for i in range(5)], "log2FoldChange"] > 1.5).all()
    assert (result.loc[[f"g{i}" for i in range(5)], "padj"] < .01).all()
    assert normalized.shape == counts.shape

