"""
Convert epidemic forecasts into JSON files for online viewing.
"""

import datetime
import errno
import h5py
import json
import logging
import numpy as np
import os.path
import pypfilt.config
import pypfilt.sweep

from . import declarative


fs_fmt = '%Y-%m-%d %H:%M:%S'
d_fmt = '%Y-%m-%d'


def dtime(dstr):
    """Convert datetime strings to datetime instances."""
    if isinstance(dstr, bytes):
        dstr = dstr.decode()
    return datetime.datetime.strptime(dstr, fs_fmt)


def date_str(dstr):
    """Convert datetime strings to date strings."""
    if isinstance(dstr, bytes):
        dstr = dstr.decode()
    dt = datetime.datetime.strptime(dstr, fs_fmt).date()
    return dt.strftime(d_fmt)


def update_obs_model(f, hdf5_file, om_dict):
    """
    Record the observation model parameters, and check that they're consistent
    across all of the input files.

    :param f: The file object from which to read the simulation output.
    :param hdf5_file: The corresponding file name for ``f``.
    :param om_dict: A dictionary of observation model parameter names/values.

    :raises ValueError: if more than one observation model is found. Note that
        the observation model may have different parameter **values** across
        all of the files.
    """
    logger = logging.getLogger(__name__)
    obs_models = list(f['meta']['param']['obs'].values())
    n_obs_models = len(obs_models)
    if n_obs_models != 1:
        raise ValueError("Found {} observation models".format(n_obs_models))
    om = obs_models[0]
    if om_dict:
        # An observation model has already been recorded, check that the
        # observation model in this file is consistent with it.
        for key in om.keys():
            if key not in om_dict:
                # A new observation model parameter has appeared.
                logger.warning("New parameter {} in {}".format(
                    key, os.path.basename(hdf5_file)))
                continue
            ok = om_dict[key] == om[key][()].item()
            if not ok:
                logger.warning("Param {} differs".format(key))
    else:
        # Record the observation model parameters.
        for key in om.keys():
            om_dict[key] = om[key][()].item()


def most_recent_obs_date(f):
    """
    Return the time of the most recent observation (a ``datetime`` instance).

    :param f: The file object from which to read the simulation output.

    :raises ValueError: if more than one observation model is found.
    """
    obs_units = list(f['data']['obs'].keys())
    if len(obs_units) != 1:
        raise ValueError("Found {} observation models".format(
            len(obs_units)))
    obs = f['data']['obs'][obs_units[0]][()]
    return max(dtime(row['date']) for row in obs)


def update_forecasts(f, hdf5_file, fs_dict, ci):
    """
    Record the forecast credible intervals and return the forecasting dates
    for which other simulation outputs should be calculated.

    :param f: The file object from which to read the simulation output.
    :param hdf5_file: The corresponding file name for ``f``.
    :param fs_dict: A dictionary of forecast credible intervals.
    :param ci: The credible intervals to record.
    """
    logger = logging.getLogger(__name__)
    # Extract the forecast credible intervals.
    fs = f['data']['forecasts'][()]
    conds = tuple(fs['prob'] == p for p in ci)
    keep = np.logical_or.reduce(conds)
    fs = fs[keep]
    # Note that forecast dates are date-time strings (%Y-%m-%d %H:%M:%S).
    fs_dates = np.unique(fs['fs_date'])
    # Check that this table contains the desired credible intervals.
    ci_levels = np.unique(fs['prob'])
    if len(ci_levels) < len(ci):
        msg = "expected CIs: {}; only found: {}"
        expect = ", ".join(str(p) for p in sorted(ci))
        found = ", ".join(str(p) for p in sorted(ci_levels))
        logger.warning(msg.format(expect, found))
    # Ignore the estimation run, if present.
    sim_end = max(dtime(dstr) for dstr in fs['date']).strftime(fs_fmt)
    if len(fs_dates) == 1 and fs_dates[0] == sim_end:
        # If the file only contains the result of an estimation run,
        # inform the user and keep these results --- they can result
        # from directly using pypfilt.run() to produce forecasts.
        last_obs = most_recent_obs_date(f)
        logger.warning('Estimation run, set fs_date = {} for {}'.format(
            last_obs.strftime(fs_fmt), os.path.basename(hdf5_file)))
        # Replace fs_date with the date of the most recent observation.
        last_obs_bs = last_obs.strftime(fs_fmt).encode()
        fs_dates = [last_obs_bs]
        fs['fs_date'] = last_obs_bs
        # Discard all rows prior to the (effective) forecasting date.
        dates = np.array([dtime(row['date']) for row in fs])
        fs = fs[dates >= last_obs]
        # Note: these files may contain duplicate rows.
        # So identify the first duplicate row (if any) and crop.
        for (n, rix) in enumerate(np.where(fs['date'] == last_obs_bs)[0]):
            # If the nth row for the date on which the forecast begins isn't
            # the nth row of the entire table, it represents the start of the
            # duplicate data, so discard all subsequent rows.
            if n != rix:
                fs = fs[:rix]
                break
    else:
        fs_dates = [d for d in fs_dates if d != sim_end]
    # Store the forecast credible intervals.
    ci_levels = np.unique(fs['prob'])
    for fs_date in fs_dates:
        mask = fs['fs_date'] == fs_date
        if not isinstance(mask, np.ndarray) or mask.shape[0] != fs.shape[0]:
            raise ValueError('Invalid fs_date comparison; {} == {}'.format(
                type(fs['fs_date'][0]), type(fs_date)))
        fs_rows = fs[mask]
        ci_dict = {}
        for ci in ci_levels:
            ci_rows = fs_rows[fs_rows['prob'] == ci]
            ci_dict[str(ci)] = [
                {"date": date_str(date),
                 "ymin": ymin,
                 "ymax": ymax}
                for (_, _, _, date, _, ymin, ymax) in ci_rows]
        fs_dict[date_str(fs_date)] = ci_dict
    # Return the forecast date(s) that should be considered for this file.
    return fs_dates


def update_peak_timing(f, hdf5_file, pkt_dict, ci, fs_dates):
    """
    Record the peak timing credible intervals.

    :param f: The file object from which to read the simulation output.
    :param hdf5_file: The corresponding file name for ``f``.
    :param pkt_dict: A dictionary of peak timing credible intervals.
    :param ci: The credible intervals to record.
    :param fs_dates: The forecasting dates for which the observations should
        be recorded.

    :raises ValueError: if more than one observation model is found. Note that
        the observation model may have different parameter **values** across
        all of the files.
    """
    logger = logging.getLogger(__name__)
    # Extract the peak timing credible intervals.
    try:
        pk = f['data']['peak_cints'][()]
    except KeyError:
        # If this table is not present, return an empty array with the
        # minimal set of required columns.
        logger.warning("No 'peak_cints' table: {}".format(
            os.path.basename(hdf5_file)))
        return
    conds = tuple(pk['prob'] == p for p in ci)
    keep = np.logical_or.reduce(conds)
    pk = pk[keep]
    ci_levels = np.unique(pk['prob'])
    if len(ci_levels) < len(ci):
        msg = "expected CIs: {}; only found: {}"
        expect = ", ".join(str(p) for p in sorted(ci))
        found = ", ".join(str(p) for p in sorted(ci_levels))
        logger.warning(msg.format(expect, found))
    for fs_date in fs_dates:
        mask = pk['fs_date'] == fs_date
        if not isinstance(mask, np.ndarray) or mask.shape[0] != pk.shape[0]:
            raise ValueError('Invalid fs_date comparison; {} == {}'.format(
                type(pk['fs_date'][0]), type(fs_date)))
        pk_rows = pk[mask]
        ci_dict = {}
        for ci in ci_levels:
            ci_rows = pk_rows[pk_rows['prob'] == ci]
            ci_dict[str(ci)] = [
                {"date": date_str(fs_date),
                 "ymin": date_str(tmin),
                 "ymax": date_str(tmax)}
                for (_, _, _, _, _smin, _smax, tmin, tmax) in ci_rows]
        pkt_dict[date_str(fs_date)] = ci_dict


def update_obs(f, hdf5_file, obs_dict, fs_dates):
    """
    Record the observations provided at each of the forecasting dates.

    :param f: The file object from which to read the simulation output.
    :param hdf5_file: The corresponding file name for ``f``.
    :param obs_dict: A dictionary of observations.
    :param fs_dates: The forecasting dates for which the observations should
        be recorded.

    :raises ValueError: if more than one observation model is found. Note that
        if the parameter **values** differ across the input files, a warning
        message will be printed but the files will still be processed.
    """
    obs_units = list(f['data']['obs'].keys())
    n_obs_units = len(obs_units)
    if n_obs_units != 1:
        raise ValueError("Found {} observation models".format(n_obs_units))
    obs = f['data']['obs'][obs_units[0]][()]
    cols = obs.dtype.names
    bs_cols = [c for c in cols if obs.dtype[c].kind == 'S']
    for fs_date in fs_dates:
        # NOTE: variable-length string values do not have an .item() method,
        # so we need to check whether the .item() method exists.
        obs_list = [
            {c: obs_row[c].item()
             if hasattr(obs_row[c], 'item') else obs_row[c]
             for c in cols}
            for obs_row in obs]
        for o in obs_list:
            # Convert byte string to Unicode strings.
            for c in bs_cols:
                o[c] = o[c].decode()
            # Ensure the date is stored as 'YYYY-MM-DD'.
            o['date'] = date_str(o['date'])
        obs_dict[date_str(fs_date)] = obs_list


def to_json(data_files, scenario, out_file, replace, plot_settings, ci=None):
    """
    Convert a set of epidemic forecasts into a JSON file for online viewing.

    :param declarative.Config config: The forecast configuration.
    :param Sequence[str] data_files: A list of forecast files (HDF5).
    :param str scenario: The forecasting scenario identifier.
    :param str out_file: The output file name.
    :param bool replace: Whether to replace (overwrite) an existing JSON file,
        rather than updating it with the provided forecasts.
    :param declarative.ForecastPlot plot_settings: The plot settings.
    :param Sequence[int] ci: The credible intervals to record (default:
        ``[0, 50, 95]``).
    """
    logger = logging.getLogger(__name__)

    if ci is None:
        ci = [0, 50, 95]

    json_data = {
        'obs': {},
        'forecasts': {},
        'timing': {},
        'obs_model': {},
        # NOTE: these are the only settings that need to be retrieved from the
        # configuration ...
        'location': scenario,
        'location_name': plot_settings['scenario_label'],
        'obs_axis_lbl': plot_settings['axis_label'],
        'obs_axis_prec': plot_settings['axis_precision'],
        'obs_datum_lbl': plot_settings['datum_label'],
        'obs_datum_prec': plot_settings['datum_precision'],
    }

    if (not replace) and os.path.isfile(out_file):
        try:
            with open(out_file, encoding='utf-8') as f:
                json_data = json.load(f)
        except json.JSONDecodeError:
            logger.warning('Could not read file "{}"'.format(out_file))
            logger.warning('Will overwrite this with a new file')

    # NOTE: sort file names to ensure the output is deterministic.
    for fname in sorted(data_files):
        with h5py.File(fname, 'r') as f:
            update_obs_model(f, fname, json_data['obs_model'])
            fs_dates = update_forecasts(f, fname, json_data['forecasts'], ci)
            update_peak_timing(f, fname, json_data['timing'], ci, fs_dates)
            update_obs(f, fname, json_data['obs'], fs_dates)

    out_dir = os.path.dirname(out_file)
    if out_dir and not os.path.isdir(out_dir):
        # Create with mode -rwxr-x---.
        try:
            logger.info('Creating {}'.format(out_dir))
            os.makedirs(out_dir, mode=0o750)
        except OSError as e:
            # Potential race condition with multiple script instances.
            if e.errno != errno.EEXIST:
                logger.warning('Could not create {}'.format(out_dir))
                raise

    logger.debug("Writing {}".format(out_file))
    with open(out_file, encoding='utf-8', mode='w') as f:
        json.dump(json_data, f, ensure_ascii=False,
                  sort_keys=True, indent=None, separators=(',', ':'))


def parser():
    """Return the command-line argument parser for ``epifx-decl-json``."""
    parser = declarative.common_parser(scenarios=False, config=True)

    op = parser.add_argument_group('Output arguments')
    op.add_argument(
        '-o', '--output', action='store', type=str, default='output.json',
        help='The name of the JSON output file')
    op.add_argument(
        '-r', '--replace', action='store_true',
        help='Replace the output file (default: update if it exists)')

    rp = parser.add_argument_group('Required arguments')
    rp.add_argument(
        '-s', '--scenario', action='store', type=str, default=None,
        help='The scenario to which the forecasts pertain')
    rp.add_argument(
        'files', metavar='HDF5_FILE', type=str, nargs='*',
        help='Forecast data file(s)')

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

    logger.info('Reading configuration from {}'.format(args['config']))
    configs = [pypfilt.config.from_file(filename)
               for filename in args['config']]

    if args['scenario'] is None:
        p.print_help()
        return 2

    scenario = args['scenario']
    forecasts = [fs for fs in pypfilt.sweep.forecasts(configs)
                 if fs.scenario_id == scenario]

    if len(forecasts) != 1:
        msg = "ERROR: invalid scenario: {}".format(args['scenario'])
        print(msg)
        return 2

    obs_models = forecasts[0].observation_models
    if len(obs_models) != 1:
        raise ValueError('Can only plot for 1 observation model')
    plot_settings = list(obs_models.values())[0].plot_settings

    to_json(data_files=args['files'], scenario=scenario,
            out_file=args['output'], replace=args['replace'],
            plot_settings=plot_settings)
