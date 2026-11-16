"""Read and validate the exact sample/feature layout of a saved MOFA model.

This module performs input validation only. It never fits a model, chooses new
features, imputes measurements or writes into the public repository.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

VIEWS = ("wes", "rnaseq", "proteomics")
QC_FILES = {
    "wes": "wes_binary_qc.csv",
    "rnaseq": "rnaseq_delta_zscore.csv",
    "proteomics": "proteomics_log2_zscore.csv",
}


def labels(values):
    return [value.decode() if isinstance(value, bytes) else str(value) for value in values]


def require_axis(observed, expected, name, allow_extra=False):
    """Fail on ambiguous/missing names rather than align arrays by position."""
    observed, expected = labels(observed), labels(expected)
    if not expected or len(expected) != len(set(expected)):
        raise ValueError(f"Empty or duplicate expected {name}")
    if len(observed) != len(set(observed)):
        raise ValueError(f"Duplicate observed {name}")
    missing = set(expected) - set(observed)
    extra = set(observed) - set(expected)
    if missing or (extra and not allow_extra):
        raise ValueError(f"{name}: missing={len(missing)}, unexpected={len(extra)}")
    return expected


def align_matrix(matrix, features, samples, view):
    """Keep only saved features, in saved order, without restandardisation."""
    import numpy as np

    matrix = matrix.copy()
    matrix.index = labels(matrix.index)
    matrix.columns = labels(matrix.columns)
    rows = require_axis(matrix.index, features, f"{view} features", allow_extra=True)
    columns = require_axis(matrix.columns, samples, f"{view} samples")
    aligned = matrix.loc[rows, columns].astype(float)
    values = aligned.to_numpy()
    if not np.isfinite(values).all():
        raise ValueError(f"Non-finite values in {view}; no silent imputation is allowed")
    if view == "wes" and not np.isin(values, [0, 1]).all():
        raise ValueError("WES training input must remain binary, not centred")
    return aligned


def read_model_layout(model_path):
    import h5py
    import numpy as np
    import pandas as pd

    with h5py.File(model_path, "r") as model:
        samples = labels(model["samples/all"][()])
        features = {view: labels(model[f"features/{view}"][()]) for view in VIEWS}
        scores = np.asarray(model["expectations/Z/all"][()]).T
    if scores.ndim != 2 or scores.shape[0] != len(samples):
        raise ValueError("Saved factor-score shape does not match stored sample names")
    require_axis(samples, samples, "model samples")
    for view, names in features.items():
        require_axis(names, names, f"model {view} features")
    factors = pd.DataFrame(scores, index=samples,
                           columns=[f"factor_{i + 1}" for i in range(scores.shape[1])])
    if not np.isfinite(scores).all():
        raise ValueError("Saved factor scores contain non-finite values")
    return samples, features, factors


def load_saved_inputs(run_root):
    import pandas as pd

    run_root = Path(run_root)
    samples, features, factors = read_model_layout(run_root / "mofa/mofa_model.hdf5")
    matrices = {}
    for view in VIEWS:
        matrix = pd.read_csv(run_root / "qc" / QC_FILES[view], index_col=0)
        matrices[view] = align_matrix(matrix, features[view], samples, view)
    return matrices, factors


def write_input_identity(destination, matrices, reference_factors):
    """Record private input identities alongside a future stability run."""
    record = {
        "reference": "saved primary MOFA model",
        "sample_order": list(reference_factors.index),
        "factor_names": list(reference_factors.columns),
        "views": {},
    }
    for view, frame in matrices.items():
        record["views"][view] = {
            "features": list(frame.index),
            "n_samples": frame.shape[1],
            "n_features": frame.shape[0],
            "matrix_sha256": hashlib.sha256(frame.to_csv().encode()).hexdigest(),
        }
    # Exclusive creation protects an existing run's provenance.
    with Path(destination).open("x", encoding="utf-8") as handle:
        json.dump(record, handle, indent=2)
        handle.write("\n")
