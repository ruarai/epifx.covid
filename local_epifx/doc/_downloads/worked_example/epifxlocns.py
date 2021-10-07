"""Local settings for epidemic forecasting."""

import epifx
import pypfilt

import epifx.summary
import pypfilt.summary

import datetime
import errno
import os.path


def get_start_date(year):
    """Begin forecasts on the first of January."""
    return datetime.datetime(year, 1, 1)


def get_until_date(year):
    """End forecasts on the last day of December."""
    return datetime.datetime(year, 12, 31)


def get_live_fs_dates(year, all_obs, fs_from):
    """Generate a forecast at the time of each observation."""
    fs_dates = sorted(o['date'] for o in all_obs)
    if fs_from is None:
        # Return only the most recent forecasting date.
        # Important: *must* be a list.
        return [fs_dates[-1]]
    else:
        return [d for d in fs_dates if d >= fs_from]


def get_scan_fs_dates(year, all_obs):
    """Generate a forecast at the time of each observation."""
    return sorted(o['date'] for o in all_obs)


def make_summary(params, all_obs, **kwargs):
    """Return a summary object for the forecast simulations."""
    meta_builder = epifx.summary.Metadata()
    meta = meta_builder.build(params)

    # Create a summary object with no tables, can now choose what tables
    # to add. Alternatively, can call epifx.summary.make().
    summary = pypfilt.summary.HDF5(params, all_obs, meta=meta, **kwargs)

    peaks = epifx.summary.PeakMonitor()
    # The summary tables to add.
    tables = [
        pypfilt.summary.ModelCIs(),
        epifx.summary.ObsLikelihood(),
        # The following tables are needed for the interactive forecast plots.
        epifx.summary.PredictiveCIs(peaks),
        epifx.summary.PeakForecastCIs(peaks),
    ]
    # Also add a summary table for each set of observations.
    # These tables are needed for the interactive forecast plots.
    for obs_unit in params['obs']:
        tables.append(epifx.summary.Obs(obs_unit))

    summary.add_tables(*tables)

    return summary


def get_locn_params(locn_settings):
    """Return the parameters dictionary for a forecast location."""
    model = epifx.SEIR()
    time = pypfilt.Datetime()
    popn_size = locn_settings['popn']
    # The number of particles.
    px_count = 2000
    # The seed for the pseudo-random number generator.
    prng_seed = 3001

    params = epifx.default_params(px_count, model, time, popn_size, prng_seed)

    # Customise the default parameters as desired.
    params['resample']['threshold'] = 0.25
    params['resample']['regularisation'] = True
    params['epifx']['stoch'] = False

    # Restrict R0 to values between 1.3 and 1.5.
    params['param_min'][4] = 1.3
    params['param_max'][4] = 1.5
    # Keep eta fixed at 1 (i.e., enforce homogeneous mixing).
    params['param_min'][7] = 1.0
    params['param_max'][7] = 1.0
    # Introduce the first infection after 12-16 weeks.
    params['param_min'][9] = 84.0
    params['param_max'][9] = 112.0
    # Update the model priors.
    params['prior'] = params['model'].priors(params)

    # Set the output and temporary directories.
    out_dir = locn_settings['out_dir']
    tmp_dir = locn_settings['tmp_dir']
    for reqd_dir in [out_dir, tmp_dir]:
        if not os.path.isdir(reqd_dir):
            # Create the directory (and missing parents) with mode -rwxr-x---.
            try:
                os.makedirs(reqd_dir, mode=0o750)
            except OSError as e:
                # Potential race condition with multiple script instances.
                if e.errno != errno.EEXIST:
                    print("Warning: could not create {}".format(reqd_dir))
                    print(e)

    params['out_dir'] = out_dir
    params['tmp_dir'] = tmp_dir

    return params


def local_settings(locn_id=None):
    """Return the forecast settings for a particular location."""
    # This dictionary will contain all of the forecast settings.
    locations = {
        'some-city': {
            'name': 'Some City',
            'popn': 1000000,
            'out_dir': '.',
            'tmp_dir': './tmp',
            'get_params': get_locn_params,
            'obs_model': epifx.obs.PopnCounts('Weekly Cases', 7),
            'obs_file': 'case-counts.ssv',
            'from_file_args': {
                'date_col': 'date',
                'value_col': 'cases',
            },
            'forecast': {
                'bg_obs': {2016: [6, 8], 2017: 8},
                'bg_var': 6,
                'pr_obs': [0.01, 0.02],
                'disp': 100,
            },
            'scan': {
                'bg_obs': [6, 8],
                'bg_var': 6,
                'pr_obs': [0.005, 0.010, 0.015, 0.020, 0.025],
                'disp': [10, 100, 1000],
            },
            'scan_years': [2015, 2016],
            'om_format': {
                'bg_obs': '02.0f',
                'bg_var': '02.0f',
                'pr_obs': '0.2f',
                'disp': '03.0f',
            },
            'om_name': {
                'bg_obs': 'bg',
                'bg_var': 'bgvar',
                'pr_obs': 'pr',
                'disp': 'disp',
            },
            'extra_args': {
                'start': get_start_date,
                'until': get_until_date,
                'live_fs_dates': get_live_fs_dates,
                'scan_fs_dates': get_scan_fs_dates,
                'make_summary': make_summary,
            },
            'json_dir': './www/data',
            'obs_axis_lbl': 'Weekly Cases',
            'obs_axis_prec': 0,
            'obs_datum_lbl': 'cases/week',
            'obs_datum_prec': 0,
        },
    }

    # Ensure each location includes its own unique ID.
    for locn in locations:
        locations[locn]['id'] = locn

    if locn_id is None:
        # Return a list of forecast locations.
        return(locations.keys())
    elif locn_id in locations:
        # Return the settings for this location.
        return(locations[locn_id])
    else:
        raise ValueError("Invalid location '{}'".format(locn_id))
