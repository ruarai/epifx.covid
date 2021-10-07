"""Test cases for the pypfilt.resample method."""

import numpy as np
import pytest

from pypfilt.resample import resample
from pypfilt.context import Scaffold


@pytest.mark.parametrize("method", ['basic', 'stratified', 'deterministic'])
def test_resample(method):
    params = {
        'resample': {
            'method': method,
            'regularisation': False,
        },
    }
    component = {
        'random': {
            'resample': np.random.RandomState(),
        },
    }

    ctx = Scaffold(params=params, component=component)

    weights = np.array([0.50, 0.25, 0.1, 0.1, 0.02, 0.02, 0.01])
    ixs = np.zeros(weights.shape)
    n_tries = 10
    for i in range(n_tries):
        x = np.array([weights, ixs]).T
        resample(ctx, x)
        if any(x[:, 1] == 0) and any(x[:, 1] == 1) and any(x[:, 1] > 1):
            # The first particle (weight 0.5) has been selected at least once,
            # as has the second particle (weight 0.25) and at least one of the
            # other particles (net weight 0.25).
            break
    else:
        pytest.fail("Failure after {} loops with {} resampling".format(
            n_tries, method))
