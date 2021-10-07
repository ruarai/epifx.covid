"""Select particles according to desired target distributions."""

import abc
import epifx.summary
import logging
import numpy as np
import pypfilt.summary
import scipy.stats


class Target(abc.ABC):
    """The base class for target particle distributions."""

    @abc.abstractmethod
    def define_summary_components(self, params):
        """
        Add summary monitors and tables so that required summary statistics
        are recorded for each proposed particle.

        :param params: The simulation parameters.
        """
        pass

    @abc.abstractmethod
    def logpdf(self, ctx, output):
        """
        Return the log of the target probability density for each particle.

        :param ctx: The simulation context.
        :param output: The state object returned by ``pypfilt.pfilter.run``;
            summary tables are located at ``output['summary'][table_name]``.
        """
        pass


class TargetPeakMVN(Target):
    """A multivariate normal distribution for the peak timing and size."""

    def __init__(self, peak_sizes, peak_times):
        """
        :param peak_sizes: An array of previously-observed peak sizes.
        :param peak_time: An array of previously-observed peak times.
        """
        exp_size = np.mean(peak_sizes)
        std_size = np.std(peak_sizes, ddof=1)
        exp_time = np.mean(peak_times)
        std_time = np.std(peak_times, ddof=1)
        self.pdf_size = scipy.stats.norm(loc=exp_size, scale=std_size)
        self.pdf_time = scipy.stats.norm(loc=exp_time, scale=std_time)
        self.log_p_max = (self.pdf_size.logpdf(exp_size) +
                          self.pdf_time.logpdf(exp_time))

    def define_summary_components(self, params):
        exp_obs_mon = pypfilt.summary.ExpectedObsMonitor()
        peak_mon = epifx.summary.PeakMonitor('exp_obs')

        params['component']['summary_monitor'] = {
            'exp_obs': exp_obs_mon,
            'peak': peak_mon,
        }
        params['component']['summary_table'] = {
            'peak_ensemble':
            epifx.summary.PeakForecastEnsembles('peak', fs_only=False),
        }

    def logpdf(self, ctx, output):
        logger = logging.getLogger(__name__)
        t = ctx.component['time']
        tbl = output['summary']['peak_ensemble']
        size = tbl['value']
        time = np.array([t.to_scalar(t.from_dtype(bs)) for bs in tbl['date']])
        logger.debug('Peak sizes: {} to {}'.format(
            np.min(size), np.max(size)))
        logger.debug('Peak times: {} to {}'.format(
            np.min(time), np.max(time)))
        log_p_size = self.pdf_size.logpdf(size)
        log_p_time = self.pdf_time.logpdf(time)
        return log_p_size + log_p_time - self.log_p_max


class TargetAny(Target):
    """A distribution that accepts all proposals with equal likelihood."""

    def define_summary_components(self, params):
        exp_obs_mon = pypfilt.summary.ExpectedObsMonitor()
        peak_mon = epifx.summary.PeakMonitor('exp_obs')

        params['component']['summary_monitor'] = {
            'exp_obs': exp_obs_mon,
            'peak': peak_mon,
        }
        params['component']['summary_table'] = {
            'peak_ensemble':
            epifx.summary.PeakForecastEnsembles('peak', fs_only=False),
        }

    def logpdf(self, ctx, output):
        return np.zeros(ctx.params['hist']['px_count'])


class Proposal(abc.ABC):
    """The base class for proposal particle distributions."""

    @abc.abstractmethod
    def sample(self, params, hist, prng):
        """
        Draw particle samples from the proposal distribution.

        :param params: The simulation parameters.
        :param hist: The particle history matrix into which the samples should
            be written.
        :param prng: The PRNG instance to use for any random sampling.
        """
        pass


class DefaultProposal(Proposal):
    """
    A proposal distribution that independently samples each parameter from the
    prior distributions provided in the simulation parameters.
    """

    def sample(self, ctx, hist, prng):
        prev_prng = ctx.component['random']['model']
        ctx.component['random']['model'] = prng
        ctx.component['model'].init(ctx, hist)
        ctx.component['random']['model'] = prev_prng


def select(params, proposal, target, seed):
    """
    Select particles according to a target distribution.

    :param params: The simulation parameters (note: the parameter dictionary
        will be emptied once the particles have been selected).
    :param proposal: The proposal distribution.
    :param target: The target distribution.
    :param seed: The PRNG seed used for sampling and accepting particles.

    :returns: The initial state vector for each accepted particle.
    :rtype: numpy.ndarray
    """
    logger = logging.getLogger(__name__)

    # Define the required summary monitors and tables.
    target.define_summary_components(params)

    # Identify the simulation period.
    start = params['time']['start']
    until = params['time']['until']

    # Create the simulation context.
    ctx = pypfilt.context.Context(params)

    px_count = params['hist']['px_count']
    prng = np.random.default_rng(seed)
    saved_hist = None
    saved_ix = 0

    while True:
        # Initialise the summary object.
        ctx.component['summary'].initialise(ctx)
        # Run the estimation pass.
        state = pypfilt.pfilter.run(ctx, start, until, [])

        # Create the history matrix for accepted particles.
        if saved_hist is None:
            saved_hist = np.zeros(state['hist'].shape[1:])
            saved_ix = 0

        # Decide which of the proposed samples to accept.
        log_pr = target.logpdf(ctx, state)
        thresh = prng.uniform(size=px_count)
        accept = np.log(thresh) < log_pr

        # Log the number of proposed particles that were accepted.
        msg = "Accept {:5d} of {:5d}"
        logger.debug(msg.format(np.sum(accept), px_count))

        # Record the accepted samples.
        upper_ix = saved_ix + np.sum(accept)
        if upper_ix > px_count:
            # Too many accepted samples, only retain a subset.
            upper_ix = px_count
            num_to_accept = upper_ix - saved_ix
            accept = np.logical_and(accept,
                                    np.cumsum(accept) <= num_to_accept)

        saved_hist[saved_ix:upper_ix, :] = state['hist'][0, accept, :]
        saved_ix = upper_ix

        if saved_ix >= px_count:
            break

    # Return the initial state vector of each accepted particle.
    state_cols = ctx.component['model'].state_size()
    return saved_hist[:, :state_cols]
