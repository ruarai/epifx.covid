"""Particle filter core: simulate time-steps and adjust particle weights."""

import logging
import numpy as np

from . import cache, check, resample
from . import obs as obs_mod
from . import state as state_mod


def reweight(ctx, hist, hist_ix, obs, max_back=None):
    """Adjust particle weights in response to some observation(s).

    :param params: The simulation parameters.
    :param hist: The particle history matrix.
    :param hist_ix: The index of the current time-step in the history matrix.
    :param obs: The observation(s) that have been made.
    :param max_back: The number of time-steps into the past when the most
        recent resampling occurred (i.e., how far back the current particle
        ordering is guaranteed to persist; default is ``None``, no limit).

    :returns: A tuple; the first element (*bool*) indicates whether resampling
        is required, the second element (*float*) is the **effective** number
        of particles (i.e., accounting for weights).
    """
    # Calculate the log-likelihood of obtaining the given observation, for
    # each particle.
    logs = obs_mod.log_llhd_of(ctx, hist, hist_ix, obs, max_back)

    # Scale the log-likelihoods so that the maximum is 0 (i.e., has a
    # likelihood of 1) to increase the chance of smaller likelihoods
    # being within the range of double-precision floating-point.
    logs = logs - np.max(logs)
    # Calculate the effective number of particles, prior to reweighting.
    prev_eff = 1.0 / sum(w * w for w in hist[hist_ix, :, -2])
    # Update the current weights.
    hist[hist_ix, :, -2] *= np.exp(logs)
    ws_sum = np.sum(sorted(hist[hist_ix, :, -2]))
    # Renormalise the weights.
    hist[hist_ix, :, -2] /= ws_sum
    if np.any(np.isnan(hist[hist_ix, :, -2])):
        # Either the new weights were all zero, or every new non-zero weight
        # is associated with a particle whose previous weight was zero.
        nans = np.sum(np.isnan(hist[hist_ix, :, -2]))
        raise ValueError("{} NaN weights; ws_sum = {}".format(nans, ws_sum))
    # Determine whether resampling is required.
    num_eff = 1.0 / sum(w * w for w in hist[hist_ix, :, -2])
    req_resample = (num_eff / ctx.params['size']
                    < ctx.params['resample']['threshold'])

    # Detect when the effective number of particles has greatly decreased.
    eff_decr = num_eff / prev_eff
    if (eff_decr < 0.1):
        # Note: this could be mitigated by replacing the weights with their
        # square roots (for example) until the decrease is sufficiently small.
        logger = logging.getLogger(__name__)
        logger.debug("Effective particles decreased by {}".format(eff_decr))

    return (req_resample, num_eff)


def __log_step(ctx, when, do_resample, num_eff=None):
    """Log the state of the particle filter when an observation is made or
    when particles have been resampled.

    :param when: The current simulation time.
    :param do_resample: Whether particles were resampled at this time-step.
    :type do_resample: bool
    :param num_eff: The effective number of particles (default is ``None``).
    :type num_eff: float
    """
    logger = logging.getLogger(__name__)
    resp = {True: 'Y', False: 'N'}
    if num_eff is not None:
        logger.debug('{} RS: {}, #px: {:7.1f}'.format(
            ctx.component['time'].to_unicode(when), resp[do_resample],
            num_eff))
    elif do_resample:
        logger.debug('{} RS: {}'.format(
            ctx.component['time'].to_unicode(when), resp[do_resample]))


def step(ctx, hist, hist_ix, step_num, when, step_obs, max_back, is_fs):
    """Perform a single time-step for every particle.

    :param params: The simulation parameters.
    :param hist: The particle history matrix.
    :param hist_ix: The index of the current time-step in the history matrix.
    :param step_num: The time-step number.
    :param when: The current simulation time.
    :param step_obs: The list of observations for this time-step.
    :param max_back: The number of time-steps into the past when the most
        recent resampling occurred; must be either a positive integer or
        ``None`` (no limit).
    :param is_fs: Indicate whether this is a forecasting simulation (i.e., no
        observations).
        For deterministic models it is useful to add some random noise when
        estimating, to allow identical particles to differ in their behaviour,
        but this is not desirable when forecasting.

    :return: ``True`` if resampling was performed, otherwise ``False``.
    """
    # Ensure we have received the entire particle history matrix.
    check.is_entire_matrix(ctx.params, hist)

    d_t = ctx.params['dt']

    # Allocate an array that enumerates the particles, if it isn't present.
    if ctx.params['px_range'] is None:
        raise ValueError('Parameter px_range not defined')

    # Matrices of previous and current state vectors.
    sc = ctx.params['hist']['state_cols']
    prev = hist[hist_ix - 1, :, 0:sc]
    curr = hist[hist_ix, :, 0:sc]

    # Step each particle forward by one time-step.
    ctx.component['model'].update(ctx, when, d_t, is_fs, prev, curr)

    # Copy the particle weights from the previous time-step.
    # These will be updated by ``reweight`` as necessary.
    hist[hist_ix, :, -2] = hist[hist_ix - 1, :, -2]

    # The particle ordering is (as yet) unchanged.
    # This will be updated by ``resample`` as necessary.
    hist[hist_ix, :, -1] = ctx.params['px_range']

    # Account for observations, if any.
    num_eff = None
    do_resample = False
    if step_obs:
        do_resample, num_eff = reweight(ctx, hist, hist_ix, step_obs,
                                        max_back)

    __log_step(ctx, when, do_resample, num_eff)

    # Perform resampling when required.
    if do_resample:
        curr = hist[hist_ix, :, 0:sc]
        ctx.component['model'].pre_resample(ctx, when, curr)


        resample.resample(ctx, hist[hist_ix])
        __log_step(ctx, when, True, ctx.params['size'])
    # Indicate whether resampling occurred at this time-step.
    return do_resample


def run(ctx, start, end, streams, state=None,
        save_when=None, save_to=None):
    """Run the particle filter against any number of data streams.

    :param ctx: The simulation parameters.
    :type ctx: pypfilt.context.Context
    :param start: The start of the simulation period.
    :param end: The (**exclusive**) end of the simulation period.
    :param streams: A list of observation streams.
    :param state: A previous simulation state as returned by, e.g., this
        function.
    :param save_when: Times at which to save the particle history matrix.
    :param save_to: The filename for saving the particle history matrix.

    :returns: The resulting simulation state: a dictionary that contains the
        simulation parameters (``'params'``), the particle history matrix
        (``'hist'``), and the summary statistics (``'summary'``).
    """
    sim_time = ctx.component['time']
    sim_time.set_period(start, end, ctx.params['steps_per_unit'])
    steps = sim_time.with_observations(*streams)
    # Determine whether this is a forecasting run, by checking whether there
    # are any observation streams.
    is_fs = not streams
    # We allow the history matrix to be provided in order to allow, e.g., for
    # forecasting from any point in a completed simulation.
    if state is None:
        hist = state_mod.history_matrix(ctx, sim_time)
        offset = 0
    else:
        hist = state['hist']
        offset = state['offset']
        # Ensure that the number of particles is recorded as a parameter.
        if hist.shape[1] != ctx.params['size']:
            raise ValueError('Invalid history matrix size')
        # Ensure the history matrix structure is appropriate.
        state_size = ctx.component['model'].state_size()
        extra = hist.shape[-1] - state_size
        if extra < 2:
            raise ValueError("Too few extra columns: {} < 2".format(extra))
        else:
            if ctx.params['hist']['state_cols'] != state_size:
                raise ValueError('Invalid state columns')
            if ctx.params['hist']['extra_cols'] != extra:
                raise ValueError('Invalid extra columns')
            # Ensure we have received the entire particle history matrix.
            check.is_entire_matrix(ctx.params, hist)

    # Allocate space for the summary statistics.
    summary = ctx.component['summary']
    summary.allocate(ctx, start, end, forecasting=is_fs)

    # Define key time-step loop variables.
    # The start of the next interval that should be summarised.
    win_start = start
    # The time of the previous time-step (if any).
    most_recent = None
    # The time-step number of the most recent resampling (if any).
    last_rs = None
    # The index of the current time-step in the state matrix.
    hist_ix = None

    # Simulate each time-step.
    for (step_num, when, obs) in steps:
        hist_ix = step_num + offset
        # Check whether the end of the history matrix has been reached.
        # If so, shift the sliding window forward in time.
        if hist_ix == hist.shape[0]:
            # Calculate summary statistics in blocks.
            if most_recent is not None:
                # The current simulation has covered a well-defined block of
                # the history matrix.
                summary.summarise(hist, sim_time, win_start, most_recent,
                                  offset)
            else:
                # If most_recent is None, no time-steps have been simulated.
                # This means, e.g., a forecasting simulation has begun at the
                # final time-step in the matrix; the correct response is to
                # calculate summary statistics only for this one time-step.
                summary.summarise(hist, sim_time, win_start, win_start,
                                  offset)
            # NOTE: win_start is the start of the next interval that will be
            # summarised. Since the current time-step is not evaluated until
            # after the above call to summary.summarise(), the next summary
            # window should start at this time-step.
            win_start = when
            shift = (ctx.params['hist']['wind_shift']
                     * ctx.params['steps_per_unit'])
            offset -= shift
            hist_ix = step_num + offset
            # Shift the sliding window forward.
            hist[:-shift, :, :] = hist[shift:, :, :]
            hist[hist_ix, :, :-2] = 0

        # Determine how many time-steps back the most recent resampling was.
        if last_rs is None:
            max_back = None
        else:
            max_back = (step_num - last_rs)

        # Simulate the current time-step.
        resampled = step(ctx, hist, hist_ix, step_num, when, obs,
                         max_back, is_fs)

        # Check whether to save the particle history matrix to disk.
        # NOTE: the summary object may not have summarised the model state
        # recently, or even at all if we haven't reached the end of the
        # sliding window at least once.
        if save_when is not None and save_to is not None:
            if when in save_when:
                # First, summarise up to the previous time-step.
                # NOTE: we do not want to summarise the current time-step,
                # because simulations that resume from this saved state will
                # begin summarising from their initial state, which is this
                # current time-step. So if the summary window starts at the
                # current time-step, we should not summarise before saving.
                if win_start < when and most_recent is None:
                    summary.summarise(hist, sim_time, win_start, win_start,
                                      offset)
                elif win_start < when:
                    summary.summarise(hist, sim_time, win_start, most_recent,
                                      offset)

                # Update the start of the next summary interval.
                win_start = when

                # Note: we only need to save the current matrix block!
                cache.save_state(save_to, ctx, when,
                                 offset=np.int32(hist_ix),
                                 hist=np.float64(hist))

        # Finally, update loop variables.
        most_recent = when
        if resampled:
            last_rs = step_num

    if hist_ix is None:
        # There were no time-steps.
        return None

    # Calculate summary statistics for the remaining time-steps.
    if most_recent is not None:
        summary.summarise(hist, sim_time, win_start, most_recent, offset)

    # Return the complete simulation state.
    return {'params': ctx.params, 'hist': hist,
            'offset': hist_ix,
            'summary': summary.get_stats()}
