"""Manage particle state vectors and their history."""

import logging
import numpy as np


def history_matrix(ctx, sim_time):
    """
    Allocate a particle history matrix of sufficient size to store an entire
    particle filter simulation.

    :param params: The simulation parameters.
    :type params: dict
    :param sim_time: The simulation period.
    :type sim_time: :py:class:`~pypfilt.Time`

    :returns: A particle history matrix.
    :rtype: numpy.ndarray
    """
    # Ensure sufficient columns to record particle weights and parents.
    state_cols = ctx.params['hist']['state_cols']
    extra_cols = ctx.params['hist']['extra_cols']
    if extra_cols < 2:
        raise ValueError("Too few extra columns: {} < 2".format(extra_cols))
    num_col_fns = len(ctx.params['hist']['extra_col_fns'])
    if num_col_fns != (extra_cols - 2):
        raise ValueError("Expected {} column functions, found: {}".format(
            extra_cols - 2, num_col_fns))
    # Determine the number of particles and their initial weights.
    px_count = ctx.params['hist']['px_count']
    # Ensure there is a strictly-positive number of particles.
    if px_count < 1:
        raise ValueError("Too few particles: {}".format(px_count))
    init_weight = 1.0 / px_count
    # Record the number of particles.
    logger = logging.getLogger(__name__)
    logger.debug("Size = {}".format(px_count))
    # Determine the number of time-steps for which to allocate space.
    wind_size = ctx.params['hist']['wind_size'] > 0
    wind_shift = ctx.params['hist']['wind_shift'] > 0
    if wind_size and wind_shift:
        num_steps = (ctx.params['hist']['wind_size']
                     * ctx.params['steps_per_unit'] + 1)
    else:
        num_steps = sim_time.step_count() + 1
    # Allocate the particle history matrix and record the initial states.
    hist = np.zeros((num_steps, px_count, state_cols + extra_cols))
    logger.debug("Hist.nbytes = {}".format(hist.nbytes))
    ctx.component['model'].init(ctx, hist[0, :, 0:state_cols])
    hist[0, :, -2] = init_weight
    hist[0, :, -1] = ctx.params['px_range']
    if extra_cols > 2:
        curr_coll = state_cols
        for name, fn in ctx.params['hist']['extra_col_fns'].items():
            fn(ctx, hist[0, :, curr_coll])
            curr_coll += 1
    # Return the allocated (and initialised) particle history matrix.
    return hist


def earlier_states(hist, ix, steps):
    """
    Return the particle states at a previous time-step, ordered with respect
    to their current arrangement.

    :param hist: The particle history matrix.
    :param ix: The current time-step index.
    :param steps: The number of steps back in time.
    """
    parent_ixs = np.arange(hist.shape[1])
    # Don't go too far back (negative indices jump into the future).
    steps = min(steps, ix)
    for i in range(steps):
        parent_ixs = hist[ix - i, parent_ixs, -1].astype(int)
    return hist[ix - steps, parent_ixs, :]
