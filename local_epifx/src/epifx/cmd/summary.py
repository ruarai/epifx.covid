"""Summarise completed scans of observation models."""

import h5py
import importlib
import logging
import numpy as np
import numpy.lib.recfunctions as recfunctions
import os
import os.path

from . import json
from . import settings


def year_from_fs_date(s, fs_date):
    """Determine the calendar year associated with one or more dates."""
    logger = logging.getLogger(__name__)

    if 'extra_args' in s and all(n in s['extra_args']
                                 for n in ['start', 'until']):
        # Determine the year by comparing forecasting dates against each of
        # the defined simulation periods.
        years = np.zeros(fs_date.shape)
        min_yr = min(fs_date).year - 1
        max_yr = max(fs_date).year + 1
        for year in range(min_yr, max_yr + 1):
            dmin = s['extra_args']['start'](year)
            dmax = s['extra_args']['until'](year)
            years[np.logical_and(fs_date >= dmin, fs_date <= dmax)] = year
        if np.any(years == 0):
            msg = '{} forecasts outside of the simulation period'.format(
                np.sum(years == 0))
            logger.error(msg)
            for ix in np.where(years == 0)[0]:
                logger.error('  -> {}'.format(fs_date[ix]))
            raise ValueError('forecasts outside of the simulation period')
        return years
    else:
        # Use the calendar year associated with the forecasting date.
        return np.array([d.year for d in fs_date])


def get_time_scale(scale_name):
    try:
        scale_mod, scale_class = scale_name.rsplit('.', 1)
        ScaleClass = getattr(importlib.import_module(scale_mod), scale_class)
        return ScaleClass()
    except AttributeError:
        raise ValueError("Time scale '{}' does not exist".format(scale_name))
    except ImportError:
        raise ValueError("Time scale import error: '{}'".format(scale_mod))
    except TypeError:
        raise ValueError("Could not instantiate '{}'".format(scale_name))


def load_tables(locn_settings, data_dir):
    logger = logging.getLogger(__name__)
    data_files = settings.find_scan_files(locn_settings, data_dir)
    logger.info("Reading forecasts from {} files in {}".format(
        len(data_files), data_dir))
    tbl_names = ['peak_size_acc', 'peak_time_acc', 'peak_cints', 'obs_llhd']
    tbls = {n: [] for n in tbl_names}
    time_scales = set()
    scale = None
    for data_file in data_files:
        with h5py.File(data_file, "r") as f:
            # Instantiate the simulation time scale.
            scale_name = f['/meta/param/time'][()]
            time_scales.add(scale_name)
            n_scales = len(time_scales)
            if n_scales > 1:
                msg = 'Multiple time scales: {}'.format(time_scales)
                raise ValueError(msg)
            elif n_scales < 1:
                msg = "No time scale for '{}'".format(data_file)
                raise ValueError(msg)
            elif scale is None:
                scale = get_time_scale(scale_name)
            # Inspect the observation model(s).
            obs_models = f['/meta/param/obs']
            om_keys = obs_models.keys()
            if len(om_keys) != 1:
                raise ValueError("{} observation models".format(len(om_keys)))
            # Note: in Python 3, h5py group methods such as keys(), values(),
            # and items() return view-like objects that cannot be sliced or
            # indexed like lists, but which support iteration.
            om = next(obs_models.values().__iter__())
            om_params = {"id_{}".format(k): om[k][()] for k in om}
            om_names = sorted(k for k in om_params)
            fields = om_names
            # Must use floats, or we can encounter "invalid type promotion".
            # om_dtypes = [om_params[k].dtype for k in om_names]
            om_values = np.array([om_params[k] for k in om_names])
            ncols = len(om_names)
            for tbl in tbl_names:
                dset = f['data'][tbl][()]
                if tbl == 'obs_llhd':
                    # Reduce to a single row, containing the mean likelihood
                    # of every observation in the estimation run.
                    fs_date = np.array([scale.from_dtype(bs)
                                        for bs in dset['fs_date']])
                    min_fs = np.min(fs_date)
                    max_fs = np.max(fs_date)
                    mask = fs_date == max_fs
                    mean_llhd = np.mean(dset['llhd'][mask])
                    year = year_from_fs_date(locn_settings,
                                             np.array([min_fs]))
                    llhd_dtype = [('year', np.int),
                                  ('llhd', dset.dtype['llhd'])]
                    dset = np.array([(year[0], mean_llhd,)], dtype=llhd_dtype)
                nrows = dset.shape[0]
                # Tag with observation model parameters.
                om_data = np.multiply(om_values[:, None],
                                      np.ones((ncols, nrows)))
                dset = recfunctions.append_fields(dset, fields, om_data,
                                                  dtypes=np.float)
                tbls[tbl].append(dset)
    for n in tbl_names:
        tbl_list = tbls[n]
        tbls[n] = np.hstack(tbl_list)
    logger.info("Finished reading forecasts")
    return tbls, scale


def using_sample_counts(tbl):
    sample_count_vars = {'id_bg_obs', 'id_disp', 'id_k_obs'}
    popn_count_vars = {'id_bg_obs', 'id_disp', 'id_pr_obs'}
    tags = {n for n in tbl.dtype.names
            if n.startswith('id_')}
    if tags >= sample_count_vars:
        return True
    elif tags >= popn_count_vars:
        return False
    else:
        raise ValueError('Unknown observation model: {}'.format(tags))


def obs_peaks(obs):
    # Identify the peak for different types of observation models.
    peaks = {}
    for o in obs:
        # Allow observations to identify their year (i.e., season).
        year = o.get('year', o['date'].year)
        if year in peaks:
            if o['value'] > peaks[year]['value']:
                peaks[year] = o
        else:
            peaks[year] = o
    return peaks


def auc_peak(s, scale, tbl, peaks, toln, weeks=7, fmt='%Y-%m-%d %H:%M:%S'):
    logger = logging.getLogger(__name__)

    peak_dates = {year: peaks[year]['date'] for year in peaks}
    peak_vals = {year: peaks[year]['value'] for year in peaks}

    tbl = tbl[tbl['toln'] == toln]

    fs_date = np.array([scale.from_dtype(bs) for bs in tbl['fs_date']])
    mask = np.zeros(tbl.shape, np.bool)
    for year in peak_dates:
        pk = peak_dates[year]
        st = scale.add_scalar(pk, -7 * weeks)
        mask = np.logical_or(mask,
                             np.logical_and(fs_date >= st, fs_date <= pk))
    tbl = tbl[mask]
    fs_date = fs_date[mask]
    years = year_from_fs_date(s, fs_date)

    # The number of rows we will produce is simply the total number of rows in
    # ``tbl`` divided by (weeks + 1).
    n_rows = tbl.shape[0] // (weeks + 1)
    if (n_rows * (weeks + 1)) != tbl.shape[0]:
        raise ValueError('How many rows for {}?'.format(tbl.shape[0]))

    # Remove the 'id_' prefix from column names.
    tags = sorted([n for n in tbl.dtype.names if n.startswith('id_')])
    cols = [('year', np.int)] + [(tag[3:], np.float) for tag in tags]
    cols += [('acc', np.float)]

    aucs = np.zeros(n_rows, dtype=cols)
    row_ix = 0
    for year in peak_vals:
        year_tbl = tbl[year == years]
        if year_tbl.shape[0] == 0:
            logger.debug("{} scans: 0".format(year))
            continue
        # Iterate over each observation model.
        seen = set()
        for row in year_tbl:
            comb = []
            for ix, tag in enumerate(tags):
                cix = year_tbl.dtype.names.index(tag)
                comb.append(row[cix])
            comb = tuple(comb)
            if comb in seen:
                continue
            else:
                seen.add(comb)
            mask = np.ones(year_tbl.shape, np.bool)
            for ix, tag in enumerate(tags):
                mask = np.logical_and(mask, year_tbl[tag] == comb[ix])
            row = [year] + list(comb) + [np.mean(year_tbl['acc'][mask])]
            aucs[row_ix] = tuple(row)
            row_ix += 1
        logger.debug("{} scans: {}".format(year, len(seen)))

    if row_ix != n_rows:
        logger.error("Produced {} rows, expected {}".format(row_ix, n_rows))
    return aucs


def best_aucs(tbl, variable, output='acc'):
    logger = logging.getLogger(__name__)
    cols = [n for n in tbl.dtype.names if n != variable and n != output]
    rows = tbl.getfield(np.dtype({n: tbl.dtype.fields[n] for n in cols}))
    names = list(rows.dtype.names) + [variable, output]
    fmt_dict = dict(tbl.dtype.descr)
    fmts = [fmt_dict[n] for n in names]
    uniq_rows = list({tuple(list(row) + [0, 0]) for row in rows})
    uniq = np.core.records.fromrecords(uniq_rows, names=names, formats=fmts)
    uniq = uniq[np.lexsort(tuple(uniq[n]
                                 for n in reversed(uniq.dtype.names)))]
    # Record the value of ``variable`` that yields the maximum ``output``.
    for ix, row in enumerate(uniq):
        comb = dict(zip(cols, row))
        mask = np.ones(tbl.shape, np.bool)
        for col, val in comb.items():
            mask = np.logical_and(mask, tbl[col] == val)
        max_val = np.max(tbl[output][mask])
        best_var = tbl[variable][mask][tbl[output][mask] == max_val]
        if best_var.shape[0] > 1:
            # Multiple matching values, so pick the first and warn the user.
            logger.warning("Found {} 'best' values for {}".format(
                best_var.shape[0], variable))
            best_var = best_var[0]
        uniq[variable][ix] = best_var
        uniq[output][ix] = max_val
    # Return the observation models that gave the best AUC scores.
    return uniq


def best_llhds(tbl):
    rows = []
    for y in sorted({y for y in tbl['year']}):
        ytbl = tbl[tbl['year'] == y]
        ytbl = ytbl[ytbl['llhd'] == np.max(ytbl['llhd'])]
        rows.extend(ytbl)
    return np.array(rows, dtype=tbl.dtype)


def best_files(tbl, locn_id, locn_settings, output='acc'):
    files = []
    om_pnames = settings.om_param_names(locn_settings, 'scan')
    for year in sorted({y for y in tbl['year']}):
        year_tbl = tbl[tbl['year'] == year]
        best_row = 0
        best_val = -1
        for ix in range(year_tbl.shape[0]):
            row_val = year_tbl[output][ix]
            if row_val > best_val:
                best_row = ix
                best_val = row_val
        om_pvals = year_tbl[om_pnames][best_row]
        om_descr = settings.descr_for(locn_settings, 'scan', om_pvals)
        out_file = settings.filename_for_scan(locn_settings, year, om_descr)
        files.append((year, out_file))
    max_len = max(len(f) for (y, f) in files)
    cols = [('year', np.int), ('file', 'S{}'.format(max_len))]
    return np.array(files, dtype=cols)


def build(s):
    logger = logging.getLogger(__name__)
    locn_id = s['id']
    observer = s['obs_model']
    obs = observer.from_file(s['obs_file'], **s['from_file_args'])
    if 'obs_filter' in s:
        obs = [o for o in obs if s['obs_filter'](o)]
    data_dir = s['out_dir']
    logger.info("Data directory is {}".format(data_dir))
    peaks = obs_peaks(obs)
    tbls, time_scale = load_tables(s, data_dir)
    out_file = settings.filename_for_scan_summary(s)

    # Calculate the predictive performance for each observation model, in the
    # lead up to the observed epidemic peak.
    auc_size = auc_peak(s, time_scale, tbls['peak_size_acc'], peaks, toln=25)
    auc_time = auc_peak(s, time_scale, tbls['peak_time_acc'], peaks, toln=10)

    if using_sample_counts(tbls['peak_size_acc']):
        variable = 'k_obs'
    else:
        variable = 'pr_obs'
    # Identify the best choices of observation probability.
    best_om_size = best_aucs(auc_size, variable)
    best_om_time = best_aucs(auc_time, variable)

    # Identify the single-best observation model for each year, and the
    # name of the corresponding output file (not including the directory).
    best_size_files = best_files(best_om_size, locn_id, s)
    best_time_files = best_files(best_om_time, locn_id, s)

    # Save the true peak size and time for each year.
    dtype_date = time_scale.dtype('date')
    peak_dtype = [('year', np.int), dtype_date, ('value', np.float)]
    peak_tbl = np.array([(year, obs['date'], obs['value'])
                         for (year, obs) in peaks.items()],
                        dtype=peak_dtype)

    # Summarise the mean observation likelihoods.
    obs_llhd = tbls['obs_llhd']
    obs_llhd.dtype.names = [n[3:] if n.startswith('id_') else n
                            for n in obs_llhd.dtype.names]
    best_llhd = best_llhds(obs_llhd)
    best_llhd_files = best_files(best_llhd, locn_id, s, output='llhd')

    # Write all of these data table to the output HDF5 file.
    logger.info("Writing results to {}".format(out_file))
    with h5py.File(out_file, "w") as f:
        # Attempt to avoid tracking times (not supported by h5py < 2.2).
        kwargs = {'track_times': False}
        # Record whether tracking times have been disabled.
        try:
            f.create_dataset('hdf5_track_times', data=False, **kwargs)
        except TypeError:
            # Cannot disable tracking times (h5py < 2.2).
            logger.info("Cannot disable tracking times")
            kwargs = {}
            f.create_dataset('hdf5_track_times', data=True, **kwargs)

        # Compress and checksum the data tables.
        kwargs['compression'] = 'gzip'
        kwargs['shuffle'] = True
        kwargs['fletcher32'] = True

        f.create_dataset('true_peaks', data=peak_tbl, **kwargs)
        grp_size = f.create_group('peak_size')
        grp_size.create_dataset('auc', data=auc_size, **kwargs)
        grp_size.create_dataset('best', data=best_om_size, **kwargs)
        grp_size.create_dataset('best_files', data=best_size_files, **kwargs)
        grp_time = f.create_group('peak_time')
        grp_time.create_dataset('auc', data=auc_time, **kwargs)
        grp_time.create_dataset('best', data=best_om_time, **kwargs)
        grp_time.create_dataset('best_files', data=best_time_files, **kwargs)
        grp_time = f.create_group('obs_llhd')
        grp_time.create_dataset('llhd', data=obs_llhd, **kwargs)
        grp_time.create_dataset('best', data=best_llhd, **kwargs)
        grp_time.create_dataset('best_files', data=best_llhd_files, **kwargs)

    # Export the best peak-timing forecasts.
    if s['json_dir']:
        for (year, filename) in best_time_files:
            json_path = settings.filename_for_scan_json(s, year)
            logger.info("Converting {} to {}".format(filename, json_path))
            files = [os.path.join(data_dir, filename)]
            json.convert(files=files, most_recent=False, locn_id=locn_id,
                         out_file=json_path, replace=True, pretty=False)


def parser():
    """Return the command-line argument parser for ``epifx-summary``."""
    p = settings.common_parser()

    fg = p.add_argument_group('Forecast settings')
    jg = fg.add_mutually_exclusive_group()
    jg.add_argument(
        '-j', '--json-dir', metavar='DIR', default=None,
        help='The directory in which JSON forecasts will be written')
    jg.add_argument(
        '-n', '--no-json', action='store_true',
        help='Do not write JSON forecasts')

    return p


def main(args=None):
    """Summarise a completed scan of observation models."""
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
        return 2

    if len(args['location']) == 0:
        p.print_help()
        return 2
    else:
        invalid = [l for l in args['location'] if l not in valid_locns]
        if invalid:
            msg = "ERROR: invalid location(s): {}".format(", ".join(invalid))
            print(msg)
            return 2

    for locn_id in args['location']:
        s = settings.local(locn_id)
        if args['no_json']:
            s['json_dir'] = None
        else:
            s['json_dir'] = settings.override('json_dir', '.', args, s)
        build(s)
