import numpy as np
import tempfile


def default_params(model, time_scale, max_days, px_count, prng_seed):
    """The default particle filter parameters.

    Memory usage can reach extreme levels with a large number of particles,
    and so it may be necessary to keep only a sliding window of the entire
    particle history matrix in memory.

    :param model: The system model.
    :param time_scale: The simulation time scale.
    :param max_days: The number of contiguous days that must be kept in memory
        (e.g., the largest observation period).
    :param px_count: The number of particles.
    :param prng_seed: The seed for the pseudo-random number generators.
    """
    details = model.describe()
    p_min = [vmin for (name, smooth, vmin, vmax) in details]
    p_max = [vmax for (name, smooth, vmin, vmax) in details]
    params = {
        'resample': {
            # Resample when the effective number of particles is 25%.
            'threshold': 0.25,
            # The deterministic method is the best resampling method, see the
            # appendix of Kitagawa 1996 (DOI:10.2307/1390750).
            'method': 'deterministic',
            # Resample from the weighted discrete probability distribution,
            # rather than using a continuous approximation (regularisation).
            'regularisation': False,
            # By default, continue without regularisation if the parameter
            # covariance matrix is not positive definite.
            'regularise_or_fail': False,
            # The minimum range of values that a parameter must have in order
            # to be subject to the post-regularised particle filter.
            'reg_toln': 1e-8,
        },
        'hist': {
            # The sliding window size, in days.
            'wind_size': 2 * max_days,
            # The amount to shift the sliding window, in days.
            'wind_shift': max_days,
            # The number of particles.
            'px_count': px_count,
            # The number of extra state columns, in addition to the model
            # state vector. Note that this number must be at least 2, since
            # the matrix must store the particle weight and parent index.
            'extra_cols': 2,
            # Functions that are responsible for initialising extra state
            # columns (except for the particle weight and parent index).
            # Mapping is name -> function.
            'extra_col_fns': {},
        },
        # Use the provided PRNG seed, if any.
        'prng_seed': prng_seed,
        # Define the PRNGs that should be created.
        'random': ['resample', 'model', 'hist_extra_cols'],
        # Define the simulation time scale.
        'component': {
            'time': time_scale,
            'model': model,
            'random': {},
            'lookup': {},
            'obs': {},
            'summary_monitor': {},
            'summary_table': {},
        },
        # Simulate 5 time-steps per unit time.
        # TODO: move into params['time']
        'steps_per_unit': 5,
        # Provide only the most recent observation period (for likelihoods).
        # TODO: move into params['hist']
        'last_n_periods': 1,
        # Whether to reduce the estimation run so that it only extends to the
        # latest forecasting date.
        'minimal_estimation_run': True,
        # An array that enumerates the particles.
        'px_range': None,
        'time': {
            # The simulation period.
            'start': None,
            'until': None,
        },
        'model': {
            # The lower bounds for each model parameter.
            'param_min': np.array(p_min),
            # The upper bounds for each model parameter.
            'param_max': np.array(p_max),
            # The model prior distributions.
            'prior': {},
        },
        'data': {
            # Observations data.
            'obs': {},
            # Lookup tables.
            'lookup': {},
        },
        'summary': {
            # If ``False`` (the default) statistics are calculated from the
            # date of the first *observation*. If ``True``, statistics are
            # calculated from the very beginning of the simulation period.
            'from_first_day': False,
            # If ``False`` (the default) statistics are calculated for the
            # initial estimation simulation and for forecasting simulations.
            # If ``True``, statistics are only calculated for forecasting
            # simulations.
            'only_forecasts': False,
            'meta': {
                # Additional packages whose versions should be recorded.
                'packages': [],
            },
        },
        # Observation model parameters.
        'obs': {},
        # Event hooks.
        'hooks': {
            'log_llhd': [],
        },
        # Directory for storing output files.
        'out_dir': '.',
        # Directory for storing temporary files.
        'tmp_dir': tempfile.gettempdir(),
    }
    return params
