"""A bootstrap particle filter for epidemic forecasting."""

import datetime
import logging
import os
import os.path

from . import cache
from . import context
from . import pfilter
from . import model
from . import obs
from . import params
from . import summary
from . import sweep
from . import time
from . import version
from . import config

__package_name__ = u'pypfilt'
__author__ = u'Rob Moss'
__email__ = u'rgmoss@unimelb.edu.au'
__copyright__ = u'2014-2020, Rob Moss'
__license__ = u'BSD 3-Clause License'
__version__ = version.__version__


# Export abstract base classes from this module.
Config = config.Config
Context = context.Context
Obs = obs.Obs
Model = model.Model
Monitor = summary.Monitor
Table = summary.Table
Datetime = time.Datetime
Scalar = time.Scalar

default_params = params.default_params

# Prevent an error message if the application does not configure logging.
log = logging.getLogger(__name__).addHandler(logging.NullHandler())


def forecasts_iter(*config_files):
    """
    Iterate over the forecast settings for each scenario defined in the
    provided configuration file(s).

    :param str config_files: The forecast configuration file(s).
    :rtype: Iterator[:class:`pypfilt.sweep.Forecasts`]

    :Examples:

    >>> import pypfilt
    >>> import pypfilt.examples.predation
    >>> pypfilt.examples.predation.write_example_files()
    >>> forecast_times = [1.0, 3.0, 5.0, 7.0, 9.0]
    >>> config_file = 'predation.toml'
    >>> data_file = 'output.hdf5'
    >>> for forecast in pypfilt.forecasts_iter(config_file):
    ...     state = pypfilt.forecast(forecast.params,
    ...                              forecast.observation_streams,
    ...                              forecast_times, filename=data_file)
    """
    configs = [config.from_file(config_file) for config_file in config_files]
    return sweep.forecasts(configs)


def simulate_from_model(params, px_count=None):
    """
    Simulate observations from a model.

    :param params: The simulation parameters.
    :return: A table of simulated observations.
    :rtype: numpy.ndarray

    :Examples:

    >>> import pypfilt
    >>> import pypfilt.examples.predation
    >>> pypfilt.examples.predation.write_example_files()
    >>> config_file = 'predation.toml'
    >>> for forecast in pypfilt.forecasts_iter(config_file):
    ...     sim_obs = pypfilt.simulate_from_model(forecast.params, px_count=1)
    ...     print(sim_obs[['date', 'unit', 'value']][:4])
    [(0., 'x', 1.33489946) (0., 'y', 0.07977887) (1., 'x', 1.93397389)
     (1., 'y', 0.58160109)]
    """
    logger = logging.getLogger(__name__)

    # Ensure that simulated observations are recorded.
    # NOTE: also remove any other summary monitors and tables.
    tbl = 'simulated_obs'
    params['component']['summary_table'] = {
        tbl: summary.SimulatedObs()
    }
    params['component']['summary_monitor'] = {}
    params['component']['summary'] = summary.HDF5(params, [])

    if px_count is not None:
        original_px_count = params['hist']['px_count']
        params['hist']['px_count'] = px_count

    ctx = Context(params)
    if px_count is not None:
        params['hist']['px_count'] = original_px_count
    params = None

    start = ctx.params['time']['start']
    until = ctx.params['time']['until']
    if start is None or until is None:
        raise ValueError("Simulation period is not defined")

    # Initialise the summary object.
    ctx.component['summary'].initialise(ctx)

    logger.info("  {}  Estimating  from {} to {}".format(
        datetime.datetime.now().strftime("%H:%M:%S"), start, until))
    state = pfilter.run(ctx, start, until, [])

    return state['summary'][tbl]


def forecast(params, streams, dates, filename):
    """Generate forecasts from various dates during a simulation.

    :param params: The simulation parameters, or a simulation context.
    :type params: Union[dict, pypfilt.Context]
    :param streams: A list of observation streams.
    :param dates: The dates at which forecasts should be generated.
    :param filename: The output file to generate (can be ``None``).

    :returns: The simulation state for each forecast date.
    """

    # Ensure that there is at least one forecasting date.
    if len(dates) < 1:
        raise ValueError("No forecasting dates specified")

    if isinstance(params, Context):
        ctx = params
    else:
        ctx = Context(params)
    params = None

    start = ctx.params['time']['start']
    end = ctx.params['time']['until']
    if start is None or end is None:
        raise ValueError("Simulation period is not defined")

    # Ensure that the forecasting dates lie within the simulation period.
    invalid_fs = [ctx.component['time'].to_unicode(d) for d in dates
                  if d < start or d >= end]
    if invalid_fs:
        raise ValueError("Invalid forecasting date(s) {}".format(invalid_fs))

    logger = logging.getLogger(__name__)

    # Initialise the summary object.
    ctx.component['summary'].initialise(ctx)

    # Generate forecasts in order from earliest to latest forecasting date.
    # Note that forecasting from the start date will duplicate the estimation
    # run (below) and is therefore redundant *if* sim['end'] is None.
    forecast_dates = [d for d in sorted(dates) if d >= start]

    # Load the most recently cached simulation state that is consistent with
    # the current observations.
    sim = cache.default(ctx, forecast_dates)
    cache_file = sim['save_to']
    update = cache.load_state(cache_file, ctx, forecast_dates)
    if update is not None:
        for (key, value) in update.items():
            sim[key] = value

    # Update the forecasting dates.
    if not sim['fs_dates']:
        logger.warning("All {} forecasting dates precede cached state".format(
            len(forecast_dates)))
        return
    forecast_dates = sim['fs_dates']

    # Update the simulation period.
    if sim['start'] is not None:
        start = sim['start']
    if sim['end'] is not None:
        # Only simulate as far as the final forecasting date, then forecast.
        # Note that this behaviour may not always be desirable, so it can be
        # disabled by setting 'minimal_estimation_run' to False.
        if ctx.params['minimal_estimation_run']:
            est_end = sim['end']
    else:
        est_end = end

    # Avoid the estimation pass when possible.
    state = None
    if start < est_end:
        logger.info("  {}  Estimating  from {} to {}".format(
            datetime.datetime.now().strftime("%H:%M:%S"), start, est_end))
        state = pfilter.run(ctx, start, est_end, streams, state=sim['state'],
                            save_when=forecast_dates, save_to=sim['save_to'])
    else:
        logger.info("  {}  No estimation pass needed for {}".format(
            datetime.datetime.now().strftime("%H:%M:%S"), est_end))

    if state is None:
        # NOTE: run() may return None if est_end < (start + dt).
        state = sim
        forecasts = {}
    else:
        # Save outputs from the estimation pass.
        # NOTE: record whether this simulation resumed from a cached state.
        if sim['start'] is not None:
            state['loaded_from_cache'] = sim['start']
        forecasts = {'complete': state}

    # Ensure the dates are ordered from latest to earliest.
    for start_date in forecast_dates:
        # We can reuse the history matrix for each forecast, since all of the
        # pertinent details are recorded in the summary.
        update = cache.load_state_at_time(cache_file, ctx, start_date)
        if update is None:
            msg = 'Cache file missing entry for forecast date {}'
            raise ValueError(msg.format(start_date))
        for (key, value) in update['state'].items():
            state[key] = value

        # The forecast may not extend to the end of the simulation period.
        fs_end = end
        if 'max_forecast_ahead' in ctx.params['time']:
            max_end = ctx.component['time'].add_scalar(
                start_date, ctx.params['time']['max_forecast_ahead'])
            if max_end < fs_end:
                fs_end = max_end

        logger.info("  {}  Forecasting from {} to {}".format(
            datetime.datetime.now().strftime("%H:%M:%S"),
            ctx.component['time'].to_unicode(start_date),
            ctx.component['time'].to_unicode(fs_end)))

        fstate = pfilter.run(ctx, start_date, fs_end, [], state=state)
        fstate['loaded_from_cache'] = start_date

        forecasts[start_date] = fstate

    # Save the observations (flattened into a single list).
    obs_list = sum(streams, [])
    forecasts['obs'] = obs_list

    # Save the forecasting results to disk.
    if filename is not None:
        logger.info("  {}  Saving to:  {}".format(
            datetime.datetime.now().strftime("%H:%M:%S"), filename))
        # Save the results in the output directory.
        filepath = os.path.join(ctx.params['out_dir'], filename)
        ctx.component['summary'].save_forecasts(forecasts, filepath)

    # Remove the temporary file and directory.
    sim['clean']()

    return forecasts


def fit(params, streams, filename):
    """
    Run a single estimation pass over the entire simulation period.

    :param params: The simulation parameters, or a simulation context.
    :type params: Union[dict, pypfilt.Context]
    :param streams: A list of observation streams.
    :param filename: The output file to generate (can be ``None``).

    :returns: The simulation state for the estimation pass.

    :Examples:

    >>> import pypfilt
    >>> import pypfilt.examples.predation
    >>> pypfilt.examples.predation.write_example_files()
    >>> config_file = 'predation.toml'
    >>> data_file = 'output.hdf5'
    >>> for forecast in pypfilt.forecasts_iter(config_file):
    ...     state = pypfilt.fit(forecast.params, forecast.observation_streams,
    ...                         filename=data_file)
    """
    if isinstance(params, Context):
        ctx = params
    else:
        ctx = Context(params)
    params = None

    start = ctx.params['time']['start']
    until = ctx.params['time']['until']
    if start is None or until is None:
        raise ValueError("Simulation period is not defined")

    logger = logging.getLogger(__name__)

    # Initialise the summary object.
    ctx.component['summary'].initialise(ctx)

    logger.info("  {}  Estimating  from {} to {}".format(
        datetime.datetime.now().strftime("%H:%M:%S"), start, until))
    state = pfilter.run(ctx, start, until, streams)
    forecasts = {'complete': state}

    # Save the observations (flattened into a single list).
    obs_list = sum(streams, [])
    forecasts['obs'] = obs_list

    # Save the forecasting results to disk.
    if filename is not None:
        logger.info("  {}  Saving to:  {}".format(
            datetime.datetime.now().strftime("%H:%M:%S"), filename))
        # Save the results in the output directory.
        filepath = os.path.join(ctx.params['out_dir'], filename)
        ctx.component['summary'].save_forecasts(forecasts, filepath)

    return forecasts
