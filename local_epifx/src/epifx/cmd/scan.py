"""Perform scans of observation models over retrospective data."""

import datetime
import logging
import os
import os.path
import pypfilt

from .. import cmd
from .. import summary
from . import settings


def forecast_dates(year, obs):
    obs_dates = sorted(set(o['date'] for o in obs))

    # Pick all observation dates from mid March until the epidemic peak.
    # First, identify the time and size of the peak.
    peak_size = max(o['value'] for o in obs)
    peak_date = [o for o in obs if o['value'] == peak_size][0]['date']

    # By default, start forecasting from the third observation in March.
    # But if there are no observations in March (e.g., the surveillance system
    # may not operate all year) then start at the earliest observation date.
    march_dates = [d for d in obs_dates if d.month == 3]
    if len(march_dates) >= 3:
        dmin = [d for d in obs_dates if d.month == 3][2]
    elif len(march_dates) > 0:
        dmin = [d for d in obs_dates if d.month == 3][0]
    else:
        dmin = min(obs_dates)

    # Ensure the final forecast is no earlier than the peak, which may or
    # may not coincide with a weekly forecasting date.
    dmax = peak_date + datetime.timedelta(days=6)
    fs_dates = [d for d in obs_dates if dmin <= d <= dmax]

    return fs_dates


def run(om_params, descr, extra):
    """
    Generate forecasts for a single set of observation model parameter values.

    :param om_params: Observation model parameters.
    :param descr:  A string that describes each observation model parameter,
        for inclusion in the output file name (ignored).
    :param extra: A dictionary of extra arguments; must include ``'year'`` and
        ``'locn_settings'``.
    """
    logger = logging.getLogger(__name__)

    year = extra['year']
    locn_settings = extra['locn_settings']
    out_file = settings.filename_for_scan(locn_settings, year, descr)

    params = settings.get_params(locn_settings)

    observer = locn_settings['obs_model']
    observer.define_params(params, **om_params)
    obs_file = locn_settings['obs_file']
    obs = observer.from_file(obs_file, year=year,
                             **locn_settings['from_file_args'])
    if 'obs_filter' in locn_settings:
        obs = [o for o in obs if locn_settings['obs_filter'](o)]

    if 'start' in extra:
        start = extra['start'](year)
    else:
        start = datetime.datetime(year, 1, 1)
    if 'until' in extra:
        until = extra['until'](year)
    else:
        until = datetime.datetime(year, 12, 31)
    if 'scan_fs_dates' in extra:
        fs_dates = extra['scan_fs_dates'](year, obs)
    else:
        fs_dates = forecast_dates(year, obs)

    sim_start = datetime.datetime.now()
    logger.debug("forecast() beginning at {}".format(
        sim_start.strftime("%H:%M:%S")))

    # Allow local settings to define the summary tables.
    make = extra.get('make_summary', summary.make)
    stats = make(params, obs, first_day=True)
    pypfilt.forecast(params, start, until, [obs], fs_dates, stats,
                     filename=out_file)
    logger.debug("forecast() finishing at {}".format(
        datetime.datetime.now().strftime("%H:%M:%S")))


def scan_iter(args):
    # Increase the "niceness" of this process (Unix platforms only).
    if args['nice'] > 0:
        try:
            os.nice(args['nice'])
        except AttributeError:
            logger = logging.getLogger(__name__)
            logger.debug("Unable to adjust niceness, os.nice() missing")

    sets = args['sets']
    subset = args['subset']

    # Use zero-based counting for partitions.
    counter = 0
    subset -= 1

    # Keep track of the output files that have been produced, to ensure there
    # are no duplicates (e.g., if observation model parameters are formatted
    # with an insufficient number of digits).
    out_files = set()

    for locn_id in args['location']:
        s = settings.local(locn_id)
        s['out_dir'] = settings.override('out_dir', '.', args, s)
        s['tmp_dir'] = settings.override('tmp_dir', '.', args, s)

        # Restrict the scan to a single year, as requested.
        if args['year'] is not None:
            if args['year'] in s['scan_years']:
                years = [args['year']]
            else:
                raise ValueError("invalid year {} for '{}'".format(
                    args['year'], locn_id))
        else:
            years = s['scan_years']

        # Perform the forecasting scan.
        for year in years:
            for (om_params, descr, _) in settings.locn_scan(s, year):
                out_file = settings.filename_for_scan(s, year, descr)
                # Ensure that this output file has not yet been used.
                if out_file in out_files:
                    raise ValueError('duplicate output file {}'.format(
                        out_file))
                else:
                    out_files.add(out_file)
                # Avoid over-writing existing files, if instructed to do so.
                if args['keep_existing']:
                    out_path = os.path.join(s['out_dir'], out_file)
                    if os.path.isfile(out_path):
                        counter = (counter + 1) % sets
                        continue
                extra = {'year': year, 'locn_settings': s}
                if 'extra_args' in s:
                    extra.update(s['extra_args'])
                # Only yield simulations in the chosen subset.
                if counter == subset:
                    yield (om_params, descr, extra)
                counter = (counter + 1) % sets


def parser():
    """Return the command-line argument parser for ``epifx-scan``."""
    p = settings.common_parser()

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
    cg.add_argument(
        '--year', metavar='YEAR', type=int, default=None,
        help='Run simulations for a single year (default: all years)')

    op = p.add_argument_group(title='Output settings')
    op.add_argument(
        '--keep-existing', action='store_true',
        help='Do not overwrite existing output files')
    op.add_argument(
        '--out-dir', type=str, default=None, metavar='DIR',
        help='The directory where forecast outputs will be saved')
    op.add_argument(
        '--tmp-dir', type=str, default=None, metavar='DIR',
        help='The directory where temporary files will be saved')

    return p


def main(args=None):
    """Perform a scan of observation models over retrospective data."""
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

    if len(args['location']) == 0:
        p.print_help()
        return 2
    else:
        invalid = [l for l in args['location'] if l not in valid_locns]
        if invalid:
            msg = "ERROR: invalid location(s): {}".format(", ".join(invalid))
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

    n_proc = args['spawn']
    if n_proc < 1:
        p.error('must use at least one process')
    elif n_proc == 1:
        for details in scan_iter(args):
            run(*details)
    else:
        cmd.run_in_parallel(run, scan_iter(args), n_proc)
