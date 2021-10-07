"""Generate forecasts from live data."""

import datetime
import errno
import logging
import os
import os.path
import pypfilt

from .. import cmd
from .. import summary
from . import json
from . import settings


def forecast_dates(year, all_obs, fs_from):
    # Use yesterday rather than today as our starting point, because the data
    # for today is certainly incomplete, and so running this script on a
    # Monday should not generate a forecast for the data up to that Monday.
    today = datetime.date.today()
    if year != today.year:
        msg = 'Live forecasting only supported for {}'.format(today.year)
        raise ValueError(msg)
    today_dt = datetime.datetime(today.year, today.month, today.day)
    yesterday = today_dt - datetime.timedelta(days=1)
    last_mon = yesterday - datetime.timedelta(days=yesterday.weekday())
    last_sun = last_mon - datetime.timedelta(days=1)

    if fs_from is None:
        # Default to the most recent Sunday.
        fs_dates = [last_sun]
    else:
        # May need to perform multiple simulations.
        if fs_from.strftime('%A') != last_sun.strftime('%A'):
            # Start simulating from the first Sunday after this date.
            to_next_sun = 7 - fs_from.isoweekday()
            fs_from = fs_from + datetime.timedelta(days=to_next_sun)
        n_weeks = 1 + (last_sun - fs_from).days // 7
        fs_dates = [fs_from + datetime.timedelta(days=7*i)
                    for i in range(n_weeks)]

    return fs_dates


def run(om_params, descr, extra):
    """
    Generate forecasts for a single set of observation model parameter values.
    Run forecasts for a single location.

    :param om_params: Observation model parameters.
    :param descr:  A string that describes each observation model parameter,
        for inclusion in the output file name (ignored).
    :param extra: A dictionary of extra arguments; must include ``'year'``,
        ``'locn_settings'``, ``'from'``, ``'inc_dates'``, and
        ``'upper_bound'``.
    """
    logger = logging.getLogger(__name__)
    s = extra['locn_settings']
    year = extra['year']
    fs_from = extra['from']
    inc_dates = extra['inc_dates']
    upper_bounds = extra['upper_bound']
    observer = s['obs_model']
    obs_file = s['obs_file']

    params = settings.get_params(s)
    cache_file = settings.filename_for_cache(s, year, descr)
    params['hist']['cache_file'] = cache_file

    # Check whether an existing cache file should be removed.
    cache_path = os.path.join(params['out_dir'], cache_file)
    if 'fresh_cache' in params and params['fresh_cache']:
        try:
            logger.info('Removing cache file {} ...'.format(cache_file))
            os.remove(cache_path)
        except OSError as e:
            # Don't raise an error if the file doesn't exist.
            if e.errno != errno.ENOENT:
                raise

    if 'start' in extra:
        start = extra['start'](year)
    else:
        start = datetime.datetime(year, 3, 1)
    if 'until' in extra:
        until = extra['until'](year)
    else:
        until = datetime.datetime(year, 10, 31)

    # Load all of the observations in this simulation period.
    all_obs = [o for o in observer.from_file(obs_file, year=year,
                                             **s['from_file_args'])
               if o['date'] > start and o['date'] <= until]
    if 'obs_filter' in s:
        all_obs = [o for o in all_obs if s['obs_filter'](o)]

    if 'live_fs_dates' in extra:
        fs_dates = extra['live_fs_dates'](year, all_obs, fs_from)
    else:
        fs_dates = forecast_dates(year, all_obs, fs_from)

    # Ignore invalid forecasting dates.
    fs_dates = [d for d in fs_dates if d > start and d < until]
    if not fs_dates:
        logger.warning("No valid forecasting dates")
        return

    observer.define_params(params, **om_params)

    # Record the forecast files that are produced.
    fs_files = []

    # Generate the forecasts.
    logger.info("Forecasting from {} ...".format(fs_dates[0]))
    for fs_date in fs_dates:
        # Extract the relevant observations.
        obs = [o for o in all_obs if o['date'] <= fs_date]
        obs_dates = [o['date'] for o in obs]

        # Mark any incomplete (underestimated) observations.
        if upper_bounds is not None:
            upper = list(upper_bounds)
        else:
            upper = None
        for o in obs:
            if o['date'] in inc_dates:
                o['incomplete'] = True
                if upper is not None:
                    # Define the upper bound for this observation.
                    o['upper_bound'] = upper.pop(0)
                    m = "incomplete observation at {} with upper bound {}"
                    logger.debug(m.format(o['date'], o['upper_bound']))
                else:
                    m = "incomplete observation at {}"
                    logger.debug(m.format(o['date']))
        # Ensure that all other incomplete observations lie in the future,
        # since "--from YYYY-MM-DD" leads to multiple forecasting dates.
        # There may be fewer upper bounds than incomplete observations.
        num_inc_future = len([o for o in inc_dates if o > fs_date])
        if upper and len(upper) > num_inc_future:
            msg = "Too many upper bounds were provided: {}"
            raise ValueError(msg.format(upper_bounds))

        if fs_date not in obs_dates:
            # Warn if there is no observation for the forecasting date.
            logger.warning("No observation for forecast date {}".format(
                fs_date))
            # Don't generate the forecast.
            continue

        out_file = settings.filename_for_forecast(s, fs_date, descr)
        fs_files.append(os.path.join(params['out_dir'], out_file))

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

    # Check whether the cache file should be removed.
    if 'remove_cache' in params and params['remove_cache']:
        try:
            logger.info('Removing cache file {} ...'.format(cache_file))
            os.remove(cache_path)
        except OSError as e:
            # Don't raise an error if the file doesn't exist.
            if e.errno != errno.ENOENT:
                raise

    if s['json_dir']:
        # Update this set of online forecasts.
        out_file = settings.filename_for_forecast_json(s, year, descr)
        logger.info("Updating {}".format(out_file))
        # Update the existing JSON file, *if* it exists.
        if os.path.isfile(out_file):
            json.convert(files=fs_files, most_recent=False, locn_id=s['id'],
                         out_file=out_file, replace=False, pretty=False)
        else:
            files = settings.find_forecast_files(s, descr, params['out_dir'])
            if files:
                json.convert(files=files, most_recent=False, locn_id=s['id'],
                             out_file=out_file, replace=True, pretty=False)
            else:
                msg = "No output files found for {}".format(out_file)
                raise ValueError(msg)


def forecast_iter(args):
    # Increase the "niceness" of this process (Unix platforms only).
    if args['nice'] > 0:
        try:
            os.nice(args['nice'])
        except AttributeError:
            logger = logging.getLogger(__name__)
            logger.debug("Unable to adjust niceness, os.nice() missing")

    sets = args['sets']
    subset = args['subset']
    year = args['year']
    incomplete = args['incomplete']

    # Use zero-based counting for partitions.
    counter = 0
    subset -= 1

    for locn_id in args['location']:
        s = settings.local(locn_id)
        if args['no_json']:
            s['json_dir'] = None
        else:
            s['json_dir'] = settings.override('json_dir', '.', args, s)
        s['out_dir'] = settings.override('out_dir', '.', args, s)
        s['tmp_dir'] = settings.override('tmp_dir', '.', args, s)

        if incomplete:
            observer = s['obs_model']
            all_obs = observer.from_file(s['obs_file'], year=year,
                                         **s['from_file_args'])
            # Ignore observation dates after the most recent Sunday.
            max_date = forecast_dates(year, all_obs, None)[0]
            all_obs = [o for o in all_obs if o['date'] <= max_date]
            if incomplete < 1 or incomplete > len(all_obs):
                raise ValueError("Invalid number of incomplete observations")
            inc_obs = all_obs[-incomplete:]
            inc_dates = [o['date'] for o in inc_obs]
        else:
            inc_dates = []

        # Perform the forecasting scan.
        for (om_params, descr, _) in settings.locn_forecasts(s, year):
            extra = {'year': year, 'locn_settings': s,
                     'from': args['from'], 'inc_dates': inc_dates,
                     'upper_bound': args['upper_bound']}
            if 'extra_args' in s:
                extra.update(s['extra_args'])
            # Only yield simulations in the chosen subset.
            if counter == subset:
                yield (om_params, descr, extra)
            counter = (counter + 1) % sets


def parser():
    """Return the command-line argument parser for ``epifx-forecast``."""
    p = settings.common_parser()

    fg = p.add_argument_group('Forecast settings')
    fg.add_argument(
        '-a', '--all', action='store_true',
        help='Generate forecasts for all known locations')
    fg.add_argument(
        '-f', '--from',
        help='First forecasting date (YYYY-MM-DD)')
    fg.add_argument(
        '-y', '--year', type=int, default=datetime.date.today().year,
        help='The forecasting year (default: current year)')
    fg.add_argument(
        '-i', '--incomplete', metavar='N', type=int, default=None,
        help='The most recent N observations are underestimates')
    fg.add_argument(
        '-u', '--upper-bound', metavar='N', type=int, action='append',
        help='The upper bound, in order from oldest to newest observation')

    cg = p.add_argument_group('Computation settings')
    cg.add_argument(
        '--spawn', metavar='N', type=int, default=1,
        help='Spawn N separate processes')
    cg.add_argument(
        '--nice', type=int, default=5,
        help='Increase the process "niceness" (default: 5)')
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

    return p


def main(args=None):
    """Generate one or more forecasts from live data."""
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

    if args['all']:
        args['location'] = valid_locns

    if len(args['location']) == 0:
        p.print_help()
        return 2
    else:
        invalid = [l for l in args['location'] if l not in valid_locns]
        if invalid:
            msg = "ERROR: invalid location(s): {}".format(", ".join(invalid))
            print(msg)
            return 2

    if args['upper_bound'] and args['incomplete']:
        if len(args['upper_bound']) != args['incomplete']:
            p.error("Number of bounds does not match --incomplete")
    elif args['upper_bound']:
        p.error("Must use both --incomplete with --upper-bound")

    if args['from'] is not None:
        try:
            args['from'] = datetime.datetime.strptime(args['from'],
                                                      '%Y-%m-%d')
        except ValueError:
            p.error("Invalid forecast date '{}'".format(args['from']))

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

    n_proc = args['spawn']
    if n_proc < 1:
        p.error('must use at least one process')
    elif n_proc == 1:
        for details in forecast_iter(args):
            run(*details)
    else:
        cmd.run_in_parallel(run, forecast_iter(args), n_proc)
