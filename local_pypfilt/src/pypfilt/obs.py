"""Observation models: expected values and log likelihoods."""

import abc
import numpy as np

from . import check, state


def expect(ctx, time, unit, period, prev, curr):
    """
    Return the expected observation value :math:`\\mathbb{E}[y_t]` for every
    every particle :math:`x_t`, at one or more times :math:`t`.

    :param ctx: The simulation context.
    :param time: The simulation time(s).
    :param unit: The observation units (see classes in :py:module`data`).
    :type unit: str
    :param period: The duration of the observation period (in days).
    :type period: int
    :param prev: The state vectors at the start of the observation period(s).
    :type prev: numpy.ndarray
    :param curr: The state vectors at the end of the observation period(s).
    :type curr: numpy.ndarray
    """
    if unit in ctx.params['obs']:
        op = ctx.params['obs'][unit]
        obs_model = ctx.component['obs'][unit]
        return obs_model.expect(ctx, op, time, period, prev, curr)
    else:
        raise ValueError("Unknown observation type '{}'".format(unit))


def log_llhd_of(ctx, hist, hist_ix, obs, max_back=None):
    """Return the log-likelihood of obtaining observations from each particle.

    :param ctx: The simulation context.
    :param hist: The particle history matrix.
    :param hist_ix: The index of the current time-step in the history matrix.
    :param obs: The observation(s) that have been made.
    :param max_back: The number of time-steps into the past when the most
        recent resampling occurred (i.e., how far back the current particle
        ordering is guaranteed to persist; default is ``None``, no limit).

    :returns: An array containing the log-likelihood for each particle.
    """
    # Ensure we have received the entire particle history matrix.
    check.is_entire_matrix(ctx.params, hist)
    if 'last_n_periods' in ctx.params and ctx.params['last_n_periods'] > 1:
        # The model requires the last N observation periods.
        rng_n = range(1, ctx.params['last_n_periods'] + 1)
        periods = set(o['period'] * n for o in obs for n in rng_n)
    else:
        periods = set(o['period'] for o in obs)
    steps_per_unit = ctx.params['steps_per_unit']

    # Extract the particle histories at every relevant prior step.
    # It may or may not be necessary to use earlier_states().
    def hist_for(period):
        """Return past state vectors in the appropriate order."""
        steps_back = steps_per_unit * period
        same_ixs = max_back is None or max_back >= steps_back
        if same_ixs:
            if steps_back > hist_ix:
                # If the observation period starts before the beginning of the
                # the simulation period, the initial state should be returned.
                return hist[0]
            else:
                return hist[hist_ix - steps_back]
        else:
            return state.earlier_states(hist, hist_ix, steps_back)
    period_hists = {period: hist_for(period) for period in periods}

    # Calculate the log-likelihood of obtaining the given observation, for
    # each particle.
    # NOTE: must pass any extra cols, they may be used by observation models.
    logs = log_llhd(ctx, obs, hist[hist_ix, :, 0:-2], period_hists,
                    hist[hist_ix, :, -2])

    return logs


def log_llhd(ctx, obs, curr, hist, weights):
    """
    Return the log-likelihood :math:`\\mathcal{l}(y_t \\mid x_t)` for the
    observation :math:`y_t` and every particle :math:`x_t`.

    :param ctx: The simulation context.
    :param obs: The list of observations for the current time-step.
    :param curr: The particle state vectors.
    :param hist: The particle state histories, indexed by observation period.
    """
    log_llhd = np.zeros(curr.shape[:-1])

    for o in obs:
        unit = o['unit']
        time = o['date']

        if unit not in ctx.params['obs']:
            raise ValueError("Unknown observation type '{}'".format(unit))

        op = ctx.params['obs'][unit]
        obs_model = ctx.component['obs'][unit]

        if hasattr(obs_model, 'log_llhd'):
            obs_llhd = obs_model.log_llhd(ctx, op, time, o, curr, hist)
        elif callable(obs_model):
            obs_llhd = obs_model(ctx, op, time, o, curr, hist)
        else:
            raise ValueError('No observation model for "{}"'.format(unit))

        ctx.call_hooks('log_llhd', o, obs_llhd, weights)

        log_llhd += obs_llhd

    return log_llhd


def simulate(ctx, time, unit, period, expected, rng=None):
    """
    Return a random sample of :math:`y_t` for each particle :math:`x_t`.

    :param ctx: The simulation context.
    :param time: The simulation time.
    :param unit: The observation unit.
    :param period: The duration of the observation period (in days).
    :param expected: The expected observation value for each particle
           :math:`x_t`.
    :param rng: The (optional) random number generator to use.
    """
    if unit in ctx.params['obs']:
        op = ctx.params['obs'][unit]
        obs_model = ctx.component['obs'][unit]
        return obs_model.simulate(ctx, op, time, period, expected, rng=rng)
    else:
        raise ValueError("Unknown observation type '{}'".format(unit))


class Obs(abc.ABC):
    """
    The base class of observation models, which defines the minimal set of
    methods that are required.
    """

    @abc.abstractmethod
    def log_llhd(self, ctx, op, time, obs, curr, hist):
        """
        Return the log-likelihood :math:`\\mathcal{l}(y_t \\mid x_t)` for the
        observation :math:`y_t` and every particle :math:`x_t`.

        :param ctx: The simulation context.
        :param op: The observation model parameters dictionary.
        :param time: The current simulation time, :math:`t`.
        :param obs: An observation for the current time-step, :math:`y_t`.
        :param curr: The particle state vectors, :math:`x_t`.
        :param hist: The particle state histories, indexed by observation
            period.
        """
        pass

    @abc.abstractmethod
    def expect(self, ctx, op, time, period, prev, curr):
        """
        Return the expected observation value :math:`\\mathbb{E}[y_t]` for
        every particle :math:`x_t`, at one or more times :math:`t`.

        :param ctx: The simulation context.
        :param op: The observation model parameters dictionary.
        :param time: The simulation time(s), :math:`t`.
        :param period: The duration of the observation period (in days).
        :param prev: The state vectors at the start of the observation
            period(s), :math:`x_t`.
        :param curr: The state vectors at the end of the observation
            period(s).
        """
        pass

    @abc.abstractmethod
    def quantiles(self, ctx, op, time, mu, wt, probs):
        r"""
        Return the values :math:`y_i` that satisfy:

        .. math::

           y_i = \inf\left\{ y : p_i \le
               \sum_i w_i \cdot \mathcal{L}(y_t \le y \mid x_t^i)\right\}

        :param ctx: The simulation context.
        :param op: The observation model parameters dictionary.
        :param time: The current simulation time, :math:`t`.
        :param mu: The expected case fraction for each particle,
            :math:`\mathbb{E}(y_t)`.
        :param wt: The weight associated with each particle, :math:`w_i`.
        :param probs: The probabilities :math:`p_i`, which **must** be sorted
            in **ascending order**.
        """
        pass

    @abc.abstractmethod
    def simulate(self, ctx, op, time, period, expected, rng=None):
        """
        Return a random sample of :math:`y_t` for each particle :math:`x_t`.

        :param ctx: The simulation context.
        :param op: The observation model parameters dictionary.
        :param time: The simulation time(s), :math:`t`.
        :param period: The duration of the observation period (in days).
        :param expected: The expected observation value for each particle
           :math:`x_t`.
        :param rng: The (optional) random number generator to use.
        """
        pass

    @abc.abstractmethod
    def from_file(self, filename, time_scale):
        """
        Load observations from a space-delimited text file with column headers
        defined in the first line.

        :param filename: The file to read.
        :param time_scale: The simulation time scale.

        :return: A list of observations, ordered as per the original file, and
            the underlying data table.
        :rtype: Tuple[List[dict], numpy.ndarray]
        """
        pass


def bisect_cdf(probs, cdf_fn, bisect_fn, y_lower, y_upper):
    r"""
    Use a bisection method to estimate the values :math:`y_i` that satisfy:

    .. math::

       y_i = \inf\left\{ y : p_i \le
           \sum_i w_i \cdot \mathcal{L}(y_t \le y \mid x_t^i)\right\}

    :param probs: The probabilities :math:`p_i`, which **must** be sorted in
        **ascending order**.
    :param cdf_fn: The CDF function
        :math:`f(y) = \sum_i w_i \cdot \mathcal{L}(y_t \le y \mid x_t^i)`.
    :param bisect_fn: The bisection function ``f(a, b)`` that either returns
        the midpoint of the interval :math:`[a, b]` or ``None`` if the
        search should stop (e.g., because a tolerance has been reached).
    :param y_lower: A lower bound for all :math:`y_i`.
    :param y_upper: An upper bound for all :math:`y_i`.
    """
    # NOTE: there may be definite lower/upper limits that cannot be exceeded
    # (e.g., 0 is a definite lower limit for Poisson and negative binomial
    # observation models). So we need to trust the observation models to
    # provide valid initial bounds, rather than checking them here.
    cdf_lower = cdf_fn(y_lower)

    bounds_lower = {pr: y_lower for pr in probs}
    bounds_upper = {pr: y_upper for pr in probs}
    qtls = np.zeros(probs.shape)
    for (ix, pr) in enumerate(probs):
        if cdf_lower >= pr:
            # The lower bound is the first value to meet or exceed this
            # threshold, so we've found y_i for this quantile.
            qtls[ix] = y_lower
            continue

        # Use a binary search to find y_i this quantile.
        y_lower = bounds_lower[pr]
        y_upper = bounds_upper[pr]

        y_mid = bisect_fn(y_lower, y_upper)
        while y_mid is not None:
            cdf_mid = cdf_fn(y_mid)

            # Check if this is a good lower or upper bound for any of the
            # remaining quantiles.
            for p in probs[ix+1:]:
                if cdf_mid <= p and y_mid > bounds_lower[p]:
                    bounds_lower[p] = y_mid
                if cdf_mid >= p and y_mid < bounds_upper[p]:
                    bounds_upper[p] = y_mid

            # Identify the half of the interval in which to search.
            if cdf_mid >= pr:
                y_upper = y_mid
            if cdf_mid <= pr:
                y_lower = y_mid
            y_mid = bisect_fn(y_lower, y_upper)

        # Record the value y_i for this quantile.
        qtls[ix] = y_upper

    return qtls
