"""Generate forecasts from live data snapshots."""

import datetime
import json
import logging
import os.path
import pypfilt
import warnings

from .. import cmd
from .. import summary
from . import settings
from .json import convert, d_fmt


def run(om_params, descr, extra):
    logger = logging.getLogger(__name__)

    # Extract required arguments from the extra dictionary.
    fs_dates = extra['fs_dates']
    obs_table = extra['obs_table']
    s = extra['locn_settings']
    year = extra['year']
    no_json = extra['no_json']
    json_src = extra['json_src']

    # Don't run the forecasts if it will clobber the input file.
    json_out = settings.filename_for_forecast_json(s, year, descr)
    if not no_json and os.path.abspath(json_out) == os.path.abspath(json_src):
        logger.warning("Will not overwrite input file {}".format(json_src))
        return

    # Record the forecast files that are produced.
    fs_files = []

    # Note that these are the same defaults as in epifx.cmd.forecast.
    if 'start' in extra:
        start = extra['start'](year)
    else:
        start = datetime.datetime(year, 3, 1)
    if 'until' in extra:
        until = extra['until'](year)
    else:
        until = datetime.datetime(year, 10, 31)

    for fs_date in fs_dates:
        # Extract the observations for this forecasting date.
        obs = obs_table[fs_date]
        # Define the simulation parameters.
        params = settings.get_params(s)
        out_file = settings.filename_for_forecast(s, fs_date, descr)
        fs_files.append(os.path.join(params['out_dir'], out_file))
        observer = s['obs_model']
        observer.define_params(params, **om_params)

        for o in obs:
            # Ensure each observation has valid metadata.
            # Note that incomplete observations are handled correctly.
            o['period'] = observer.period
            o['unit'] = observer.unit
            o['source'] = json_src

        # Generate the forecast for this date.
        logger.info(" {} {}".format(fs_date, out_file))
        sim_start = datetime.datetime.now()
        logger.debug("forecast() beginning at {}".format(
            sim_start.strftime("%H:%M:%S")))
        # Allow local settings to define the summary tables.
        make = extra.get('make_summary', summary.make)
        stats = make(params, obs, first_day=True, only_fs=True)
        pypfilt.forecast(params, start, until, [obs], [fs_date], stats,
                         filename=out_file)
        logger.debug("forecast() finishing at {}".format(
            datetime.datetime.now().strftime("%H:%M:%S")))

    if no_json:
        return

    # Write the forecasts to a single JSON file.
    logger.info("Writing {}".format(json_out))
    convert(fs_files, True, s['id'], json_out, replace=True, pretty=False)


def replay_iter(args):
    """
    Return a generator that performs a parameter sweep for live forecasts,
    intended for use with ``epifx.cmd.replay.run``
    """
    logger = logging.getLogger(__name__)

    sets = args['sets']
    subset = args['subset']
    json_file = args['json_file']
    no_json = args['no_json']
    locn_id = args['location']

    try:
        with open(json_file, encoding='utf-8') as f:
            json_data = json.load(f)
    except json.JSONDecodeError as e:
        logger.error("Could not read file '{}'".format(json_file))
        raise e

    if locn_id is None:
        # By default, use the same forecasting location as that found in the
        # input JSON file.
        locn_id = json_data['location']

    # Use the settings for this forecasting location, but allow values
    # provided on the command-line to override the provided values.
    s = settings.local(locn_id)
    s['out_dir'] = settings.override('out_dir', '.', args, s)
    s['tmp_dir'] = settings.override('tmp_dir', '.', args, s)
    s['json_dir'] = settings.override('json_dir', '.', args, s)

    fs_dates = sorted([datetime.datetime.strptime(d, d_fmt)
                       for d in json_data['forecasts'].keys()])

    # Read the observations that were provided at each forecasting date.
    obs_table = {}
    for fs_date in fs_dates:
        obs_table[fs_date] = json_data['obs'][str(fs_date.date())]
        for o in obs_table[fs_date]:
            o['date'] = datetime.datetime.strptime(o['date'], d_fmt)

    # Define perfect upper-bound estimates for each observation.
    if args['perfect']:
        final_fs = fs_dates[-1]
        final_obs = obs_table[final_fs]
        perfect_tbl = {o['date']: o['upper_bound']
                       if o.get('incomplete', False) else o['value']
                       for o in final_obs}
        msg_ok = "{:%Y-%m-%d} forecast {:%Y-%m-%d} value {} upper_bound {}"
        msg_ex = "{:%Y-%m-%d} forecast {:%Y-%m-%d} value {} exceeds {}"
        for fs_date in fs_dates[:-1]:
            for o in obs_table[fs_date]:
                # (Re)define the observation upper bounds.
                final_val = perfect_tbl[o['date']]
                if o['value'] < final_val:
                    o['incomplete'] = True
                    o['upper_bound'] = final_val
                    logger.debug(msg_ok.format(fs_date, o['date'], o['value'],
                                               final_val))
                elif o['value'] > final_val:
                    # The observed value is larger than the final value.
                    # This is ignored for now, although we could adjust the
                    # observation models so that they no longer assume that
                    # the upper bound is larger than the observed value (in
                    # which case, it should also be renamed!).
                    logger.debug(msg_ex.format(fs_date, o['date'], o['value'],
                                               final_val))

    # Use zero-based counting for partitions.
    counter = 0
    subset -= 1

    # Determine the year from the first forecasting date.
    year = fs_dates[0].year

    # Perform the forecasting scan.
    for (om_params, descr, _) in settings.locn_forecasts(s, year):
        # Group all of the forecasts for a single observation model into a
        # single simulation, so that they can share a common cache file
        # without any race issues arising.
        extra = {'year': year, 'locn_settings': s, 'obs_table': obs_table,
                 'fs_dates': fs_dates, 'no_json': no_json,
                 'json_src': json_file}
        if 'extra_args' in s:
            extra.update(s['extra_args'])
        # Only yield simulations in the chosen subset.
        if counter == subset:
            yield (om_params, descr, extra)
        counter = (counter + 1) % sets


def parser():
    """Return the command-line argument parser for ``epifx-replay``."""
    p = settings.common_parser(locns=False)

    cg = p.add_argument_group('Computation settings')
    cg.add_argument(
        '--spawn', metavar='N', type=int, default=1,
        help='Spawn N separate processes')
    cg.add_argument(
        '--subset', metavar='I', type=int, default=None,
        help='Run the Ith subset of the simulations')
    cg.add_argument(
        '--sets', metavar='S', type=int, default=1,
        help='Divide the simulations into S subsets.')

    og = p.add_argument_group(title='Output settings')
    og.add_argument(
        '--out-dir', type=str, default=None, metavar='DIR',
        help='The directory where forecast outputs will be saved')
    og.add_argument(
        '--tmp-dir', type=str, default=None, metavar='DIR',
        help='The directory where temporary files will be saved')

    jg = og.add_mutually_exclusive_group()
    jg.add_argument(
        '-j', '--json-dir', metavar='DIR', default=None,
        help='The directory in which JSON forecasts will be written')
    jg.add_argument(
        '-n', '--no-json', action='store_true',
        help='Do not write JSON forecasts')

    fg = p.add_argument_group(title='Forecast settings')
    fg.add_argument(
        '-p', '--perfect-upper-bounds', action='store_true', dest='perfect',
        help='Set observation upper bounds to their final values')
    fg.add_argument(
        'json_file', type=str,
        help='The JSON file containing the live data snapshots')
    fg.add_argument(
        'location', type=str, nargs='?', default=None,
        help='The forecasting location identifier (optional)')

    return p


def main(args=None):
    """Generate forecasts from live data snapshots."""
    p = parser()
    if args is None:
        args = vars(p.parse_args())
    else:
        args = vars(p.parse_args(args))
    logging.basicConfig(level=args['loglevel'])

    try:
        valid_locns = settings.locations()
    except settings.NoLocalSettings as e:
        print(e)
        return(2)

    if args['location'] is not None and args['location'] not in valid_locns:
        msg = "ERROR: invalid location: {}".format(args['location'])
        print(msg)
        return 2

    if args['sets'] == 1:
        if args['subset'] is not None and args['subset'] != 1:
            p.error('subset cannot exceed the number of sets')
        elif args['subset'] is None:
            args['subset'] = 1
    elif args['sets'] < 1:
        p.error('the number of sets must be positive')
    elif args['subset'] is None:
        p.error('the subset must be defined when there are multiple sets')
    elif args['subset'] < 1:
        p.error('the subset must be positive')
    elif args['subset'] > args['sets']:
        p.error('the subset cannot exceed the number of sets')

    warnings.simplefilter('error', category=BytesWarning)

    n_proc = args['spawn']
    if n_proc < 1:
        p.error('must use at least one process')
    elif n_proc == 1:
        for details in replay_iter(args):
            run(*details)
    else:
        cmd.run_in_parallel(run, replay_iter(args), n_proc)
