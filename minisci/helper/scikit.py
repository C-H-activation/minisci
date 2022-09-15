"""Scikit-learn helper functions."""

import joblib
import numpy as np


def create_machine_learning_barriers(features: dict, model_path: str) -> dict:
    """Create machine-learning barriers from pre-calculated features.
    Args:
    features: physicochemical features (dict)
    Returns:
    Mapping of aromatic carbon (key) to activation barrier in kcal/mol (value).
    """
    # initialize barriers
    barriers = {}

    # load machine learning model from path (e.g., /path/to/ml.mod)
    model = joblib.load(model_path)

    # loop over features
    X = []
    for _, feature_vector in features.items():
        X.append(feature_vector)
    X = np.array(X)
    # predict barriers in kcal/mol
    yp = model.predict(X)
    counter = 0
    for carbon_idx, _ in features.items():
        barriers[carbon_idx] = yp[counter]
        counter += 1
    return barriers


def create_boltzmann_weights(
    barriers: dict, in_percent=True, T=298.15, R=0.0083145
) -> dict:
    """Create Boltzmann weights from barrieres in kcal/mol."""
    # initialize weights
    weights = dict()

    norm = 0
    scaling_factor = 4.184
    for carbon_idx, barrier in barriers.items():
        exponential = np.exp(-(barrier * scaling_factor / (R * T)))
        weights[carbon_idx] = exponential
        norm += exponential

    # normalize
    for carbon_idx, weight in weights.items():
        weights[carbon_idx] = weight / norm

    if in_percent:
        for carbon_idx, weight in weights.items():
            weights[carbon_idx] = weight * 100
    return weights
