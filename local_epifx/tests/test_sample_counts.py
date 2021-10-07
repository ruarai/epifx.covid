"""Test cases for the epifx.obs.SampleCounts class."""

import numpy as np
import pytest

from epifx.obs import SampleCounts


def test_llhd_boundaries():
    patients = 1000
    om = SampleCounts(__file__, 7, patients)
    prob = 0.01
    disp = 100
    for cases in [0, 1, 500, 999, 1000]:
        log_pm = om.logpmf(x=cases, prob=prob, size=patients, disp=disp)
        assert(log_pm < 0)
    old_err_settings = np.seterr(invalid='ignore')
    for cases in [-1, 1001]:
        # The betaln function will yield an invalid floating-point value, and
        # by default a warning message will be printed (suppressed by the call
        # to np.seterr, above).
        log_pm = om.logpmf(x=cases, prob=prob, size=patients, disp=disp)
        assert(np.isnan(log_pm) or np.isinf(log_pm))
        assert(not np.isfinite(log_pm))
    np.seterr(all='raise')
    for cases in [-1]:
        with pytest.raises(FloatingPointError):
            log_pm = om.logpmf(x=cases, prob=prob, size=patients, disp=disp)
    np.seterr(**old_err_settings)
