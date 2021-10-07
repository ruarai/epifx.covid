"""Test cases for the pypfilt.stats module."""

import numpy as np
import pytest

import pypfilt.stats


def test_cov_wt_R_example():
    """
    Ensure that cov_wt() produces the same output as the R function cov.wt(),
    for the example provided in the R documentation.
    """
    w = np.array([0, 0, 0, 1, 1, 1, 1, 1, 0, 0], dtype=float)
    x = np.c_[np.arange(1, 11), np.arange(1, 11)]
    y = pypfilt.stats.cov_wt(x, w, cor=False)
    exp_cov = np.array([[2.5, 2.5], [2.5, 2.5]])
    np.testing.assert_allclose(exp_cov, y)
    y = pypfilt.stats.cov_wt(x, w, cor=True)
    exp_cor = np.array([[1, 1], [1, 1]])
    np.testing.assert_allclose(exp_cor, y)
    x[3:7, 1] = [8, 7, 6, 5]
    y = pypfilt.stats.cov_wt(x, w, cor=False)
    exp_cov = np.array([[2.5, -0.5], [-0.5, 1.7]])
    np.testing.assert_allclose(exp_cov, y)
    y = pypfilt.stats.cov_wt(x, w, cor=True)
    exp_cor = np.array([[1, -0.2425356250363329], [-0.2425356250363329, 1]])
    np.testing.assert_allclose(exp_cor, y)


def test_cov_wt_ndim():
    """Ensure that cov_wt() fails when the weights are not one-dimensional."""
    w = np.array([[1, 2, 3], [4, 5, 6]])
    x = np.c_[np.arange(1, 11), np.arange(1, 11)]
    with pytest.raises(ValueError):
        pypfilt.stats.cov_wt(x, w)


def test_cov_wt_empty():
    """Ensure that cov_wt() fails when the weights are empty."""
    w = np.array([])
    x = np.c_[np.arange(1, 11), np.arange(1, 11)]
    with pytest.raises(ValueError):
        pypfilt.stats.cov_wt(x, w)


def test_cov_wt_too_few():
    """
    Ensure that cov_wt() fails when there are fewer weights than observations.
    """
    w = np.array([0, 0, 0, 1, 1, 1, 1, 1, 0], dtype=float)
    x = np.c_[np.arange(1, 11), np.arange(1, 11)]
    with pytest.raises(ValueError):
        pypfilt.stats.cov_wt(x, w, cor=False)


def test_cov_wt_too_many():
    """
    Ensure that cov_wt() fails when there are more weights than observations.
    """
    w = np.array([0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 1], dtype=float)
    x = np.c_[np.arange(1, 11), np.arange(1, 11)]
    with pytest.raises(ValueError):
        pypfilt.stats.cov_wt(x, w, cor=False)


def test_cov_wt_sign():
    """Ensure that cov_wt() fails when a weight is negative."""
    w = np.array([0, 0, 0, 1, 1, 1, -1, 1, 0, 0], dtype=float)
    x = np.c_[np.arange(1, 11), np.arange(1, 11)]
    with pytest.raises(ValueError):
        pypfilt.stats.cov_wt(x, w, cor=False)


def test_avg_var_wt():
    """
    Ensure avg_var_wt() outputs matches that of, e.g., the GNU Scientific
    Library (GSL) function gsl_stats_wvariance, as documented here:
    `https://www.gnu.org/software/gsl/manual/html_node/Weighted-Samples.html`
    (also see the section "Weighted sample variance" of
    `https://en.wikipedia.org/wiki/Weighted_arithmetic_mean` and the
    correction factor for reliability weights).
    """
    x = np.array([3.7, 3.3, 3.5, 2.8])
    w = np.array([5, 5, 4, 1]) / 15
    exp_mean = 259 / 75
    exp_bi = np.sum(w * (x - exp_mean) ** 2)
    exp_unbi = exp_bi / (1 - np.sum(w ** 2))
    (mean1, varn_unbi) = pypfilt.stats.avg_var_wt(x, w, biased=False)
    (mean2, varn_bi) = pypfilt.stats.avg_var_wt(x, w, biased=True)
    assert np.allclose(mean1, exp_mean), 'weighted mean'
    assert np.allclose(mean2, exp_mean), 'weighted mean'
    assert np.allclose(exp_bi, varn_bi), 'weighted biased variance'
    assert np.allclose(exp_unbi, varn_unbi), 'weighted unbiased variance'


def test_qtl_wt_unif():
    """
    Ensure qtl_wt() produces the same output as the R function
    ``wtd.quantile(x, w, probs, normwt=TRUE)``, provided by the Hmisc package.
    """
    # Values are the integers 1 to 10 (inclusive) with uniform weights.
    x = np.arange(1, 11)
    probs = np.array([0.025, 0.05, 0.10, 0.25, 0.50,
                      0.75, 0.90, 0.95, 0.975])
    w = np.ones(x.shape) / len(x)
    exp_qtls = np.array([1.225, 1.450, 1.900, 3.250, 5.500,
                         7.750, 9.100, 9.550, 9.775])
    qtls = pypfilt.stats.qtl_wt(x, w, probs)
    assert len(qtls) == len(exp_qtls)
    assert np.allclose(qtls, exp_qtls)


def test_qtl_wt_nonunif():
    """
    Ensure qtl_wt() produces the same output as the R function
    ``wtd.quantile(x, w, probs, normwt=TRUE)``, provided by the Hmisc package.
    """
    # Values are the integers 1 to 10 (inclusive) with non-uniform weights.
    x = np.arange(1, 11)
    probs = np.array([0.025, 0.05, 0.10, 0.25, 0.50,
                      0.75, 0.90, 0.95, 0.975])
    w = np.array([0.05, 0.05, 0.05, 0.1, 0.25,
                  0.25, 0.1, 0.05, 0.05, 0.05])
    exp_qtls = np.array([2.45, 2.90, 3.80, 5.0, 5.50,
                         6.75, 8.20, 9.10, 9.55])
    qtls = pypfilt.stats.qtl_wt(x, w, probs)
    assert len(qtls) == len(exp_qtls)
    assert np.allclose(qtls, exp_qtls)


def test_qlt_wt_dups():
    """
    Ensure qtl_wt() produces the same output as the R function
    ``wtd.quantile(x, w, probs, normwt=TRUE)``, provided by the Hmisc package.
    """
    # Test with duplicate values.
    x = np.arange(1, 11)
    x[:3] = 1
    x[7:] = 8
    probs = np.array([0.025, 0.05, 0.10, 0.25, 0.50,
                      0.75, 0.90, 0.95, 0.975])
    w = np.array([0.05, 0.05, 0.05, 0.1, 0.25,
                  0.25, 0.1, 0.05, 0.05, 0.05])
    exp_qtls = np.array([1.675, 2.350, 3.700, 5.000, 5.500,
                         6.750, 8.000, 8.000, 8.000])
    qtls = pypfilt.stats.qtl_wt(x, w, probs)
    assert len(qtls) == len(exp_qtls)
    assert np.allclose(qtls, exp_qtls)


def test_cred_wt_1():
    """Ensure that weighted credible intervals are calculated correctly."""
    x = np.arange(1, 11)
    w = np.array([0, 0.05, 0.1, 0.1, 0.25, 0.25, 0.1, 0.1, 0.05, 0])
    cs = [0, 50, 90]
    intervals = pypfilt.stats.cred_wt(x, w, cs)
    # Ensure that all intervals are reported.
    assert all(k in cs for k in intervals)
    assert all(k in intervals for k in cs)


def test_cred_wt_2():
    """
    Ensure that weighted credible intervals are calculated correctly when
    duplicate values are present.
    """
    x = np.array([1, 1, 1, 2, 3])
    w = np.array([0.3, 0.3, 0.3, 0.05, 0.05])
    cs = [0, 50, 90]
    intervals = pypfilt.stats.cred_wt(x, w, cs)
    # Ensure that all intervals are reported.
    assert all(k in cs for k in intervals)
    assert all(k in intervals for k in cs)
