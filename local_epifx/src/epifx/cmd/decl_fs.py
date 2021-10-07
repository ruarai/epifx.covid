"""
Run forecasts from a declarative configuration.
"""

import datetime
import errno
import logging
import os
import pypfilt
import pypfilt.config
import pypfilt.sweep
import re

from pathlib import Path

from . import declarative
from . import decl_json
from .. import cmd


def get_forecast_dates(all_obs, fs_from):
    """
    Return the dates for which forecasts should be generated.

    :param all_obs: The available observations.
    :type all_obs: List[Dict[str, Any]]
    :param fs_from: The earliest potential forecasting date; if this is
        ``None`` then only the most recent forecasting date will be returned.
    :type fs_from: Optional[datetime.datetime]
    """
    newest_obs = max(obs['date'] for obs in all_obs)
    if fs_from is None:
        return [newest_obs]
    return sorted(obs['date'] for obs in all_obs
                  if obs['date'] >= fs_from)


def run(forecast, forecast_dates):
    """
    Run forecast simulations for each forecasting date.

    :param forecast: The forecast settings for a scenario.
    :type forecast: pypfilt.sweep.Forecasts
    :param forecast_dates: The dates at which to run forecasts.
    :type forecast_dates: List[datetime.datetime]
    :return: The simulation state and forecast output file for each
        forecasting date.
    :rtype: Dict[datetime.datetime, Dict]
    """
    logger = logging.getLogger(__name__)

    params = forecast.params
    start = params['time']['start']
    until = params['time']['until']

    cache_file = filename_for_description(forecast.scenario_id,
                                          forecast.file_descriptor,
                                          prefix='cache', infix=None)
    params['hist']['cache_file'] = cache_file
    cache_path = Path(params['out_dir']) / cache_file

    if params['fresh_cache']:
        try:
            logger.info('Removing cache file {} ...'.format(cache_file))
            os.remove(cache_path)
        except OSError as e:
            # Don't raise an error if the file doesn't exist.
            if e.errno != errno.ENOENT:
                raise

    # Determine the valid forecasting dates.
    forecast_dates = [d for d in forecast_dates if d > start and d < until]

    # Run each forecast in turn.
    forecast_files = []
    forecast_states = {}
    out_dir = Path(params['out_dir'])
    for fs_date in forecast_dates:
        logger.info("Forecasting from {} ...".format(fs_date))

        if not any(fs_date == o['date'] for o in forecast.all_observations):
            # Warn if there is no observation for the forecasting date.
            logger.warning("No observation for forecast date {}".format(
                fs_date))
            # Don't generate the forecast.
            continue

        # NOTE: we don't want to include the time component in filenames, so
        # we have to handle datetime.datetime values as a special case.
        if isinstance(fs_date, datetime.datetime):
            fs_str = str(fs_date.date())
        else:
            time = params['component']['time']
            fs_str = time.to_unicode(fs_date)
        # NOTE: no need to add the out_dir prefix.
        fs_file = filename_for_description(forecast.scenario_id,
                                           forecast.file_descriptor,
                                           infix=fs_str)
        forecast_files.append(out_dir / fs_file)
        sim_start = datetime.datetime.now()
        logger.debug("forecast() beginning at {}".format(
            sim_start.strftime("%H:%M:%S")))

        state = pypfilt.forecast(params,
                                 forecast.observation_streams,
                                 [fs_date],
                                 filename=fs_file)
        state['forecast_file'] = fs_file
        forecast_states[fs_date] = state
        logger.debug("forecast() finishing at {}".format(
            datetime.datetime.now().strftime("%H:%M:%S")))

    # Check whether the cache file should be removed.
    if params['remove_cache']:
        try:
            logger.info('Removing cache file {} ...'.format(cache_file))
            os.remove(cache_path)
        except OSError as e:
            # Don't raise an error if the file doesn't exist.
            if e.errno != errno.ENOENT:
                raise

    if 'json_dir' in params:
        obs_models = forecast.observation_models
        if len(obs_models) != 1:
            raise ValueError('Can only plot for 1 observation model')
        plot_settings = list(obs_models.values())[0].plot_settings
        json_dir = params['json_dir']
        json_file = filename_for_description(forecast.scenario_id,
                                             forecast.file_descriptor,
                                             ext='.json')
        json_path = Path(json_dir) / json_file
        logger.info("Updating {}".format(json_path))
        if os.path.isfile(json_path):
            # Update an existing JSON file.
            decl_json.to_json(data_files=forecast_files,
                              scenario=forecast.scenario_id,
                              out_file=json_path,
                              replace=False,
                              plot_settings=plot_settings)
        else:
            # If there is no existing JSON file, include all forecast outputs
            # for this forecasting configuration in the output.
            # NOTE: this includes forecasts that were not generated by the
            # simulations that were run above.
            regex = ('^' + forecast.scenario_id + r'-\d{4}-\d{2}-\d{2}-'
                     + forecast.file_descriptor + '\\.hdf5$')
            files = [Path(params['out_dir']) / filename
                     for filename in os.listdir(params['out_dir'])
                     if re.match(regex, filename)]
            if files:
                decl_json.to_json(data_files=files,
                                  scenario=forecast.scenario_id,
                                  out_file=json_path,
                                  replace=True,
                                  plot_settings=plot_settings)
            else:
                msg = "No output files found for {}".format(json_path)
                raise ValueError(msg)

    return forecast_states


def run_multiprocess(mp_tuple, forecast_from):
    forecast = pypfilt.sweep.get_forecasts_mp(mp_tuple)
    forecast_dates = get_forecast_dates(forecast.all_observations,
                                        forecast_from)
    return run(forecast, forecast_dates)


def filename_for_description(scenario, descr, prefix=None, infix=None,
                             suffix=None, ext=None):
    if ext is None:
        ext = '.hdf5'

    fields = []
    if prefix is not None:
        fields.append(prefix)
    fields.append(scenario)
    if infix is not None:
        fields.append(infix)
    fields.append(descr)
    if suffix is not None:
        fields.append(suffix)
    return '-'.join(fields) + ext


def parser():
    """Return the command-line argument parser for ``epifx-decl-fs``."""
    parser = declarative.common_parser(scenarios=True, config=True)

    fg = parser.add_argument_group('Forecast settings')

    fg.add_argument(
        '-f', '--from',
        help='First forecasting date (YYYY-MM-DD)')

    cg = parser.add_argument_group('Computation settings')
    cg.add_argument(
        '--spawn', metavar='N', type=int, default=1,
        help='Spawn N separate processes')
    cg.add_argument(
        '--nice', type=int, default=5,
        help='Increase the process "niceness" (default: 5)')

    # NOTE: this argument currently has no effect.
    # og = parser.add_argument_group(title='Output settings')
    # og.add_argument(
    #     '-n', '--no-json', action='store_true',
    #     help='Do not write JSON forecasts')

    return parser


def main(args=None):
    """Generate forecasts from live data."""
    p = parser()
    if args is None:
        args = vars(p.parse_args())
    else:
        args = vars(p.parse_args(args))
    logging.basicConfig(level=args['loglevel'])
    logger = logging.getLogger(__name__)

    if args['config'] is None:
        p.error('must specify at least one configuration file')

    logger.info('Reading configuration from {}'.format(args['config']))
    configs = [pypfilt.config.from_file(filename)
               for filename in args['config']]

    if args['from'] is not None:
        try:
            args['from'] = datetime.datetime.strptime(args['from'],
                                                      '%Y-%m-%d')
        except ValueError:
            p.error("Invalid forecast date '{}'".format(args['from']))

    scenario_ids = args['scenario']
    forecast_from = args['from']

    n_proc = args['spawn']
    if n_proc < 1:
        p.error('must use at least one process')
    elif n_proc == 1:
        forecasts = (fs for fs in pypfilt.sweep.forecasts(configs)
                     if (not scenario_ids)
                     or fs.scenario_id in scenario_ids)
        for forecast in forecasts:
            fs_dates = get_forecast_dates(forecast.all_observations,
                                          forecast_from)
            run(forecast, fs_dates)
    else:
        forecasts = ((mp_tuple, forecast_from)
                     for mp_tuple in pypfilt.sweep.forecasts_mp(configs)
                     if (not scenario_ids)
                     or mp_tuple[2] in scenario_ids)
        cmd.run_in_parallel(run_multiprocess, forecasts, n_proc)
