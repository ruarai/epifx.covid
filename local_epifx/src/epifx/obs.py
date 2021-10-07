"""Observation models: expected values and log likelihoods."""

import abc
import numpy as np
import pypfilt.obs
import scipy.stats

from pypfilt.io import read_table
from scipy.special import betaln, gammaln


class Obs(pypfilt.obs.Obs):
    """
    The base class of observation models, which defines the minimal set of
    methods that are required.
    """

    @abc.abstractmethod
    def llhd_in(self, op, time, mu, wt, y0, y1):
        """
        Return the weighted likelihood that :math:`y_t \\in [y_0, y_1)`:

        .. math::

           \\sum_i w_i \\cdot \\mathcal{L}(y_0 \\le y_t < y_1 \\mid x_t^i)

        :param op: The observation model parameters dictionary.
        :param time: The current simulation time, :math:`t`.
        :param mu: The expected case fraction for each particle,
            :math:`\\mathbb{E}(y_t)`.
        :param wt: The weight associated with each particle, :math:`w_i`.
        :param y0: The (inclusive) minimum fraction of cases, :math:`y_0`.
        :param y1: The (exclusive) maximum fraction of cases, :math:`y_1`.
        """
        pass


class SampleCounts(Obs):
    """
    Generic observation model for relating disease incidence to count data
    where the denominator is reported.
    """

    def __init__(self, obs_unit, obs_period, denom, upper_bound_as_obs=False,
                 k_obs_lookup=None):
        """
        :param obs_unit: A descriptive name for the data.
        :param obs_period: The observation period (in days).
        :param denom: The denominator to use when calculating likelihoods and
            quantiles in the absence of an actual observation.
        :param upper_bound_as_obs: Treat upper bounds as **point estimates**.
        :param k_obs_lookup: The name of a lookup table for the
            disease-related increase in observation rate
            (:math:`\\kappa_\\mathrm{obs}`). By default, the value in the
            parameters dictionary is used.
        """
        self.unit = obs_unit
        self.period = obs_period
        self.denom = denom
        self.upper_bound_as_obs = upper_bound_as_obs
        self.k_obs_lookup = k_obs_lookup

    @staticmethod
    def logpmf(x, prob, size, disp):
        """
        Return the log of the probability mass at :math:`x`.

        :param x: The number of cases (observed numerator :math:`x`).
        :param prob: The expected fraction of all patients that are cases.
        :param size: The number of patients (observed denominator).
        :param disp: The dispersion parameter (:math:`k`).
        """
        return (gammaln(size + 1) - gammaln(x + 1) - gammaln(size - x + 1) -
                betaln(disp * (1 - prob), disp * prob) +
                betaln(size - x + disp * (1 - prob), x + disp * prob))

    @classmethod
    def interval_pmf(cls, x0, x1, prob, size, disp, log=True):
        """
        Return the (log of the) probability mass over the interval
        :math:`(x_0, x1]`.

        :param x0: The (exclusive) minimum number of cases (:math:`x_0`).
        :param x1: The (inclusive) maximum number of cases (:math:`x_1`).
        :param prob: The expected fraction of all patients that are cases.
        :param size: The number of patients (observed denominator).
        :param disp: The dispersion parameter (:math:`k`).
        :param log: Whether to return the log of the probability mass.
        """
        total = np.zeros(np.shape(prob))
        for x in range(x0 + 1, x1 + 1):
            total += np.exp(cls.logpmf(x, prob, size, disp))
        if log:
            # Handle particles with zero mass in this interval.
            total[total <= 0.0] = np.finfo(total.dtype).tiny
            return np.log(total)
        return total

    def expect(self, ctx, op, time, period, prev, curr):
        """
        Calculate the expected observation value :math:`\\mathbb{E}[y_t]` for
        every particle :math:`x_t`.
        """
        pr_inf = ctx.component['model'].pr_inf(prev, curr)
        if self.k_obs_lookup is not None:
            col_ix = ctx.params['hist']['state_cols']
            col_ix += ctx.params['sample_lookup_tables'][self.k_obs_lookup]
            ixs = curr[..., col_ix].astype(int)
            values = (ctx.component['lookup']
                      [self.k_obs_lookup].lookup(time))
            k_obs = values[ixs]
        else:
            k_obs = op['k_obs']
        return op['bg_obs'] + pr_inf * k_obs

    def simulate(self, ctx, op, time, period, pr, rng=None):
        """
        Return a random sample for each particle.
        """
        disp = self.effective_disp(pr, op)
        # To generate a random value from the beta-binomial distribution,
        # first draw 'p' randomly from the Beta(a, b) distribution, then draw
        # 'x' (numerator) from the Bin(N, p) distribution.
        alpha = pr * disp
        beta = (1 - pr) * disp
        # NOTE: using the default sample size, and returning the sample as a
        # *proportion* with respect to that sample size.
        if rng is None:
            p = scipy.stats.beta(alpha, beta).rvs()
            return scipy.stats.binom(self.denom, p).rvs() / self.denom
        else:
            p = rng.beta(alpha, beta)
            return rng.binom(self.denom, p) / self.denom

    def log_llhd(self, params, op, time, obs, curr, hist):
        """
        Calculate the log-likelihood :math:`\\mathcal{l}(y_t \\mid x_t)` for
        the observation :math:`y_t` (``obs``) and every particle :math:`x_t`.

        If it is known (or suspected) that the observed value will increase in
        the future --- when ``obs['incomplete'] == True`` --- then the
        log-likehood :math:`\\mathcal{l}(y > y_t \\mid x_t)` is calculated
        instead (i.e., the log of the *survival function*).

        If an upper bound to this increase is also known (or estimated) ---
        when ``obs['upper_bound']`` is defined --- then the log-likelihood
        :math:`\\mathcal{l}(y_u \\ge y > y_t \\mid x_t)` is calculated
        instead.

        The upper bound can also be treated as a **point estimate** by setting
        ``upper_bound_as_obs = True`` --- then the
        log-likelihood :math:`\\mathcal{l}(y_u \\mid x_t)` is calculated.
        """
        period = obs['period']
        pr = self.expect(params, op, time, period, hist[period], curr)
        num = obs['numerator']
        denom = obs['denominator']
        disp = self.effective_disp(pr, op)

        if 'incomplete' in obs and obs['incomplete']:
            if 'upper_bound' in obs:
                # Calculate the log-likelihood over the interval from the
                # observed value to this upper bound.
                num_max = obs['upper_bound']
                if self.upper_bound_as_obs:
                    # Return the likelihood of observing the upper bound.
                    return self.logpmf(num_max, pr, denom, disp)
                return self.interval_pmf(num, num_max, pr, denom, disp)
            else:
                # Calculate the log-likelihood of observing a strictly greater
                # value than reported by this incomplete observation.
                return self.logsf(num, pr, denom, disp)
        # Calculate the log-likehood of the observed value.
        return self.logpmf(num, pr, denom, disp)

    def trapz_qtls(self, mu, wt, disp, probs, trapz_w):
        """
        Approximate the quantile function using trapezoidal integration.

        :param mu: The expected case fraction for each particle.
        :param wt: The weight associated with each particle.
        :param disp: The dispersion parameter for each particle.
        :param probs: The cumulative probabilities that define the credible
            interval boundaries (must lie in the interval :math:`(0, 1)` and
            must be sorted in *ascending order*).
        :param trapz_w: The trapezoid width; quantiles will be approximated by
            linear interpolation between trapezoid boundaries.
        """
        denom = self.denom
        if not np.array_equal(probs, np.sort(probs)):
            raise ValueError("unsorted intervals: {}".format(probs))
        cints = np.zeros(probs.shape)
        # Record point and cumulative probabilities for a subset of [0, Nt].
        csums = np.zeros((denom + 1,))
        pmass = np.zeros((denom + 1,))

        # Calculate the probability mass for a specific number of cases.
        def pm_at(x):
            return np.dot(wt, np.exp(self.logpmf(x, mu, denom, disp)))

        # Calculate cumulative probabilities at the boundaries.
        pmass[0] = pm_at(0)
        csums[0] = pmass[0]
        pmass[-1] = pm_at(denom)
        csums[-1] = 1.0
        # Calculate cumulative probabilities for each trapezoid in turn.
        x_lwr = 0
        x_upr = 0
        # Trapezoids span intervals (x_lwr, x_upr], weight accordingly.
        w_lwr = (trapz_w - 1) / 2
        w_upr = (trapz_w + 1) / 2
        for pix, pr in enumerate(probs):
            while csums[x_upr] < pr:
                x_lwr = x_upr
                x_upr += trapz_w
                if x_upr > denom:
                    x_upr = denom
                    break
                pmass[x_upr] = pm_at(x_upr)
                csums[x_upr] = (csums[x_lwr] + w_lwr * pmass[x_lwr] +
                                w_upr * pmass[x_upr])
            # Have found bounds on x for the given probability.
            if x_lwr == x_upr:
                x_est = x_lwr
            else:
                # Use linear interpolation to estimate the value of x.
                slope = (csums[x_upr] - csums[x_lwr]) / (x_upr - x_lwr)
                diff = pr - csums[x_lwr]
                x_est = np.rint(x_lwr + diff / slope).astype(int)
            cints[pix] = x_est / denom
        return cints

    def trapz_approx(self, mu, wt, disp, x0, x1, n_samp):
        """
        Approximate the probability mass over the interval :math:`[x_0, x_1)`
        using a trapezoidal integration scheme.

        :param mu: The expected case fraction for each particle.
        :param wt: The weight associated with each particle.
        :param disp: The dispersion parameter for each particle.
        :param x0: The (inclusive) minimum number of cases (:math:`x_0`).
        :param x1: The (exclusive) maximum number of cases (:math:`x_1`).
        :param n_samp: The number of samples over the :math:`[x_0, x_1)`
            interval (i.e., set to :math:`t + 1` for :math:`t` trapezoids).
        """
        x = np.linspace(x0, x1 - 1, num=n_samp, dtype=int)
        grid_x, grid_mu = np.meshgrid(x, mu, copy=False)
        # Note: we use self.denom as the denominator here; this must be
        # consistent with the denominator used to calculate x0 and x1.
        probs = np.exp(self.logpmf(grid_x, grid_mu, self.denom, disp))
        widths = 0.5 * np.ediff1d(x)
        # Note: the end-points need to be scaled by (w + 1)/2 for intervals of
        # width w (x is discrete, so they both have *unit* width).
        traps = np.r_[0.5, widths] + np.r_[widths, 0.5]
        trap_probs = probs * traps[None, :]
        trap_sums = np.dot(wt, np.sum(trap_probs, axis=1))
        return trap_sums

    def effective_disp(self, mu, op, expand_dim=False):
        """
        Return the dispersion parameter for each particle, subject to an
        optional lower bound imposed on the variance.
        """
        disp = op['disp']

        if 'bg_var' in op and op['bg_var'] > 0:
            # Ensure the variance is not smaller than the variance in the
            # background signal.
            disp = op['disp'] * np.ones(mu.shape)
            min_var = op['bg_var']
            frac_var = (mu * (1 - mu) * (1 + (self.denom - 1) / (disp + 1))
                        / self.denom)
            mask_v = frac_var < min_var
            if np.any(mask_v):
                c = mu[mask_v] * (1 - mu[mask_v]) / self.denom
                disp[mask_v] = (self.denom * c - min_var) / (min_var - c)
            else:
                return op['disp']
            if expand_dim:
                disp = disp[:, None]

        return disp

    def llhd_in(self, op, time, mu, wt, y0, y1):
        """
        Return the probability mass over the interval :math:`[y_0, y1)`.

        :param op: The observation model parameters dictionary.
        :param time: The current simulation time, :math:`t`.
        :param mu: The expected fraction of all patients that are cases.
        :param wt: The weight associated with each value of ``mu``.
        :param y0: The (inclusive) minimum fraction of cases (:math:`y_0`).
        :param y1: The (exclusive) maximum fraction of cases (:math:`y_0`).
        """
        x0 = np.ceil(y0 * self.denom).astype(int)
        x1 = np.ceil(y1 * self.denom).astype(int)
        disp = self.effective_disp(mu, op, expand_dim=True)
        diff = x1 - x0
        if diff <= 20:
            # Sufficiently small interval, calculate the mass at every x.
            x = np.linspace(x0, x1 - 1, num=diff, dtype=int)
            grid_x, grid_mu = np.meshgrid(x, mu, copy=False)
            probs = np.exp(self.logpmf(grid_x, grid_mu, self.denom, disp))
            return np.dot(wt, np.sum(probs, axis=1))
        elif diff <= 50:
            n_samp = 6
        else:
            # Maximum relative error < 0.1% on CDC national data.
            n_samp = 9
        return self.trapz_approx(mu, wt, disp, x0, x1, n_samp)

    def quantiles(self, ctx, op, time, mu, wt, probs):
        r"""
        Return the observations :math:`y_i` that satisfy:

        .. math::

           y_i = \inf\left\{ y \in \mathbb{N} : p_i \le
               \sum_i w_i \cdot \mathcal{L}(y_t \le y \mid x_t^i)\right\}

        :param op: The observation model parameters dictionary.
        :param time: The current simulation time, :math:`t`.
        :param mu: The expected case fraction for each particle,
            :math:`\mathbb{E}(y_t)`.
        :param wt: The weight associated with each particle, :math:`w_i`.
        :param probs: The probabilities :math:`p_i`, which **must** be sorted
            in **ascending order**.
        """
        disp = self.effective_disp(mu, op)
        # Determine the trapezoid width by setting an upper limit on the
        # number of trapezoids that can be used to span the interval [0, 1],
        # and defining a minimum trapezoid width.
        denom = self.denom
        max_traps = 10000
        min_width = 10
        trapz_w = max(min_width, np.ceil(denom / max_traps).astype(int))
        return self.trapz_qtls(mu, wt, disp, probs, trapz_w)

    @classmethod
    def logsf(cls, x, prob, size, disp):
        """
        Return the log of the survival function
        :math:`(1 - \\mathrm{CDF}(x))`.

        :param x: The number of cases (observed numerator :math:`x`).
        :param prob: The expected fraction of all patients that are cases.
        :param size: The number of patients (observed denominator).
        :param disp: The dispersion parameter (:math:`k`).
        """
        total = np.ones(np.size(prob))
        for x in range(x + 1):
            total -= np.exp(cls.logpmf(x, prob, size, disp))
        # Handle particles with zero mass in this interval.
        total[total == 0.0] = np.finfo(total.dtype).tiny
        return np.log(total)

    def from_file(self, filename, time_scale, year=None, time_col='to',
                  value_col='cases', denom_col='patients'):
        """
        Load count data from a space-delimited text file with column headers
        defined in the first line.

        Note that returned observation *values* represent the *fraction* of
        patients that were counted as cases, **not** the *absolute number* of
        cases.
        The number of cases and the number of patients are recorded under the
        ``'numerator'`` and ``'denominator'`` keys, respectively.

        :param filename: The file to read.
        :param year: Only returns observations for a specific year.
            The default behaviour is to return all recorded observations.
        :param time_col: The name of the observation time column.
        :param value_col: The name of the observation value column (reported
            as absolute values, **not** fractions).
        :param denom_col: The name of the observation denominator column.
        :return: A list of observations, ordered as per the original file, and
            the underlying data table.

        :raises ValueError: If a denominator or value is negative, or if the
            value exceeds the denominator.
        """
        cols = [time_scale.column(time_col), (value_col, np.int32),
                (denom_col, np.int32)]
        if year is not None:
            year_col = 'year'
            cols.insert(0, (year_col, np.int32))
        df = read_table(filename, cols)

        if year is not None:
            df = df[df[year_col] == year]

        # Perform some basic validation checks.
        if np.any(df[denom_col] < 0):
            raise ValueError("Observation denominator is negative")
        elif np.any(df[value_col] < 0):
            raise ValueError("Observed value is negative")
        elif np.any(df[value_col] > df[denom_col]):
            raise ValueError("Observed value exceeds denominator")

        # Return observations with non-zero denominators.
        nrows = df.shape[0]
        obs_list = [{'date': df[time_col][i],
                     'value': df[value_col][i] / df[denom_col][i],
                     'numerator': df[value_col][i],
                     'denominator': df[denom_col][i],
                     'unit': self.unit,
                     'period': self.period,
                     'source': str(filename)}
                    for i in range(nrows) if df[denom_col][i] > 0 and
                    df[value_col][i] > 0]
        return (obs_list, df)


class PopnCounts(Obs):
    """
    Generic observation model for relating disease incidence to count data
    where the denominator is assumed or known to be the population size.
    """

    def __init__(self, obs_unit, obs_period, upper_bound_as_obs=False,
                 pr_obs_lookup=None):
        """
        :param obs_unit: A descriptive name for the data.
        :param obs_period: The observation period (in days).
        :param upper_bound_as_obs: Treat upper bounds as **point estimates**.
        :param pr_obs_lookup: The name of a lookup table for the observation
            probability (:math:`p_\\mathrm{obs}`). By default, the value in
            the parameters dictionary is used.
        """
        self.unit = obs_unit
        self.period = obs_period
        self.upper_bound_as_obs = upper_bound_as_obs
        self.pr_obs_lookup = pr_obs_lookup

    def expect(self, ctx, op, time, period, prev, curr):
        """
        Calculate the expected observation value :math:`\\mathbb{E}[y_t]` for
        every particle :math:`x_t`.
        """
        n = ctx.component['model'].population_size()
        pr_inf = ctx.component['model'].pr_inf(prev, curr)
        if self.pr_obs_lookup is not None:
            col_ix = ctx.params['hist']['state_cols']
            col_ix += ctx.params['sample_lookup_tables'][self.pr_obs_lookup]
            ixs = curr[..., col_ix].astype(int)
            values = (ctx.component['lookup']
                      [self.pr_obs_lookup].lookup(time))
            pr_obs = values[ixs]
        else:
            pr_obs = op['pr_obs']
        return (1 - pr_inf) * op['bg_obs'] + pr_inf * pr_obs * n

    def effective_disp(self, mu, op):
        """
        Return the dispersion parameter for each particle, subject to an
        optional lower bound imposed on the variance.
        """
        nb_k = op['disp']

        # Ensure the variance is not smaller than the variance in the
        # background signal.
        if 'bg_var' in op and op['bg_var'] > 0:
            nb_k = op['disp'] * np.ones(mu.shape)
            min_var = op['bg_var']
            nb_var = mu + np.square(mu) / nb_k
            mask_v = nb_var < min_var
            if np.any(mask_v):
                nb_k[mask_v] = np.square(mu[mask_v]) / (min_var - mu[mask_v])

        return nb_k

    def simulate(self, ctx, op, time, period, mu, rng=None):
        """
        Return a random sample for each particle.
        """
        nb_k = self.effective_disp(mu, op)
        nb_pr = nb_k / (nb_k + mu)
        if rng is None:
            nbinom = scipy.stats.nbinom(nb_k, nb_pr)
            return nbinom.rvs()
        else:
            return rng.negative_binomial(nb_k, nb_pr)

    def log_llhd(self, params, op, time, obs, curr, hist):
        """
        Calculate the log-likelihood :math:`\\mathcal{l}(y_t \\mid x_t)` for
        the observation :math:`y_t` (``obs``) and every particle :math:`x_t`.

        If it is known (or suspected) that the observed value will increase in
        the future --- when ``obs['incomplete'] == True`` --- then the
        log-likehood :math:`\\mathcal{l}(y > y_t \\mid x_t)` is calculated
        instead (i.e., the log of the *survival function*).

        If an upper bound to this increase is also known (or estimated) ---
        when ``obs['upper_bound']`` is defined --- then the log-likelihood
        :math:`\\mathcal{l}(y_u \\ge y > y_t \\mid x_t)` is calculated
        instead.

        The upper bound can also be treated as a **point estimate** by setting
        ``upper_bound_as_obs = True`` --- then the
        log-likelihood :math:`\\mathcal{l}(y_u \\mid x_t)` is calculated.
        """
        period = obs['period']
        mu = self.expect(params, op, time, period, hist[period], curr)
        nb_k = self.effective_disp(mu, op)
        nb_pr = nb_k / (nb_k + mu)

        if 'pr_detect' in obs:
            # NOTE: scale the probability to account for incomplete detection.
            nb_pr = nb_pr / (nb_pr + obs['pr_detect'] * (1 - nb_pr))

        nbinom = scipy.stats.nbinom(nb_k, nb_pr)
        if 'incomplete' in obs and obs['incomplete']:
            if 'upper_bound' in obs:
                if self.upper_bound_as_obs:
                    # Return the likelihood of observing the upper bound.
                    return nbinom.logpmf(obs['upper_bound'])
                # Calculate the likelihood over the interval from the observed
                # value to this upper bound, and return its logarithm.
                cdf_u = nbinom.cdf(obs['upper_bound'])
                cdf_l = nbinom.cdf(obs['value'])
                # Handle particles with zero mass in this interval.
                probs = cdf_u - cdf_l
                probs[probs <= 0] = np.finfo(probs.dtype).tiny
                return np.log(probs)
            else:
                # Return the likelihood of observing a strictly greater value
                # than the value reported by this incomplete observation.
                return nbinom.logsf(obs['value'])
        return nbinom.logpmf(obs['value'])

    def llhd_in(self, op, time, mu, wt, y0, y1):
        """
        Return the probability mass in :math:`[y_0, y1)`.

        :param op: The observation model parameters dictionary.
        :param time: The current simulation time, :math:`t`.
        :param mu: The expected case fraction for each particle.
        :param wt: The weight associated with each particle.
        :param y0: The (inclusive) minimum fraction of cases (:math:`y_0`).
        :param y1: The (exclusive) maximum fraction of cases (:math:`y_0`).
        """
        nb_k = self.effective_disp(mu, op)
        nb_pr = nb_k / (nb_k + mu)
        nbinom = scipy.stats.nbinom(nb_k, nb_pr)
        return np.dot(wt, nbinom.cdf(y1 - 1) - nbinom.cdf(y0 - 1))

    def quantiles(self, ctx, op, time, mu, wt, probs):
        r"""
        Return the observations :math:`y_i` that satisfy:

        .. math::

           y_i = \inf\left\{ y \in \mathbb{N} : p_i \le
               \sum_i w_i \cdot \mathcal{L}(y_t \le y \mid x_t^i)\right\}

        :param op: The observation model parameters dictionary.
        :param time: The current simulation time, :math:`t`.
        :param mu: The expected case fraction for each particle,
            :math:`\mathbb{E}(y_t)`.
        :param wt: The weight associated with each particle, :math:`w_i`.
        :param probs: The probabilities :math:`p_i`, which **must** be sorted
            in **ascending order**.
        """
        nb_k = self.effective_disp(mu, op)
        nb_pr = nb_k / (nb_k + mu)
        nbinom = scipy.stats.nbinom(nb_k, nb_pr)

        def cdf(y):
            return np.dot(wt, nbinom.cdf(y))

        y_min = 0
        # Find a satisfactory upper bound for y_i.
        max_pr = np.max(probs)
        for scale in np.logspace(1, 8, base=2, num=8):
            y_max = np.max(mu) * scale
            c_max = cdf(y_max)
            if c_max >= max_pr:
                # We've found a satisfactory upper bound for y_i.
                break
        else:
            msg = "Could not find a satisfactory upper bound for p = {}"
            raise ValueError(msg.format(max_pr))

        def bisect(a, b):
            if b > a + 1:
                return np.rint((a + b) / 2).astype(int)
            else:
                return None

        return pypfilt.obs.bisect_cdf(probs, cdf, bisect, y_min, y_max)

    def from_file(self, filename, time_scale, year=None, time_col='to',
                  value_col='count', ub_col=None, pr_detect_col=None):
        """
        Load count data from a space-delimited text file with column headers
        defined in the first line.

        :param filename: The file to read.
        :param year: Only returns observations for a specific year.
            The default behaviour is to return all recorded observations.
        :param time_col: The name of the observation date column.
        :param value_col: The name of the observation value column.
        :param ub_col: The name of the estimated upper-bound column, optional.
        :param pr_detect_col: The name of the column that defines a detection
            probability that accounts "incomplete" observations due to, e.g.,
            delays in reporting.
        :return: A list of observations, ordered as per the original file, and
            the underlying data table.
        """
        cols = [time_scale.column(time_col), (value_col, np.int32)]
        if year is not None:
            year_col = 'year'
            cols.insert(0, (year_col, np.int32))
        if ub_col is not None:
            cols.append((ub_col, np.int32))
        if pr_detect_col is not None:
            cols.append((pr_detect_col, np.float_))
        df = read_table(filename, cols)

        if year is not None:
            df = df[df[year_col] == year]

        nrows = df.shape[0]
        if ub_col is None:
            obs_list = [{'date': df[time_col][i],
                         'value': df[value_col][i],
                         'unit': self.unit,
                         'period': self.period,
                         'source': str(filename)}
                        for i in range(nrows)]
        else:
            obs_list = [{'date': df[time_col][i],
                         'value': df[value_col][i],
                         'incomplete': df[ub_col][i] > df[value_col][i],
                         'upper_bound': df[ub_col][i],
                         'unit': self.unit,
                         'period': self.period,
                         'source': str(filename)}
                        for i in range(nrows)]

        if pr_detect_col is not None:
            for i in range(nrows):
                obs_list[i]['pr_detect'] = df[pr_detect_col][i]

        return (obs_list, df)
