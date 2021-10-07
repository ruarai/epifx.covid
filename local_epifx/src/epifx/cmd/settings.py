"""Command-line argument parsers used by epifx tools."""

import argparse
import datetime
import itertools
import logging
import os
import os.path
import re
import sys

from .. import version


def common_parser(locns=True):
    """A command-line argument parser with common settings."""
    if locns:
        usage_str = '%(prog)s [options] location [location ...]'
    else:
        usage_str = '%(prog)s [options]'

    p = argparse.ArgumentParser(usage=usage_str, add_help=False)

    h = p.add_argument_group("Information")
    h.add_argument(
        '-h', '--help', action='help',
        help='Show this help message')
    h.add_argument("--version", action="version",
                   version="epifx {}".format(version.__version__),
                   help="Print the version information and exit.")

    og = p.add_argument_group('Output options')
    log = og.add_mutually_exclusive_group()
    log.add_argument(
        '-d', '--debug', action='store_const', dest='loglevel',
        help='Enable debugging output',
        const=logging.DEBUG, default=logging.INFO)
    log.add_argument(
        '-q', '--quiet', action='store_const', dest='loglevel',
        help='Suppress logging output',
        const=logging.WARNING)

    if locns:
        lg = p.add_argument_group('Forecasting location')
        lg.add_argument(
            'location', nargs='*', metavar='location',
            help='Forecasting location identifier(s)')

    return p


class NoLocalSettings(Exception):
    def __init__(self, locn_id=None):
        if locn_id is None:
            msg = "Could not import {0}, have you created {0}.py?"
        else:
            msg = "Could not find '{1}' in {0}.py, is it defined?"
        Exception.__init__(self, msg.format('epifxlocns', locn_id))


def local(locn_id=None):
    """
    Load locally-defined forecasting settings.

    :param locn_id: The location identifier.
    """
    logger = logging.getLogger(__name__)
    cwd_abs = os.path.abspath('.')
    add_cwd = cwd_abs not in sys.path
    try:
        if add_cwd:
            logger.debug("Adding '.' to sys.path".format(cwd_abs))
            sys.path.insert(0, cwd_abs)
            import epifxlocns
        else:
            import epifxlocns
        if hasattr(epifxlocns, '__file__'):
            mod_dir = os.path.dirname(epifxlocns.__file__)
            rel_dir = os.path.relpath(mod_dir)
            if rel_dir != '.':
                logger.debug("Locations module found in {}".format(rel_dir))
        else:
            logger.debug("Locations module has no '__file__' attribute")
        if hasattr(epifxlocns, 'local_settings'):
            if locn_id:
                return epifxlocns.local_settings(locn_id)
            else:
                return epifxlocns.local_settings()
        else:
            raise NoLocalSettings('local_settings')
    except ImportError:
        raise NoLocalSettings()
    finally:
        if add_cwd:
            logger.debug("Removing '.' from sys.path".format(cwd_abs))
            sys.path.pop(0)


def locations():
    """Return the list of valid location identifiers."""
    return local()


def get_params(settings):
    return settings['get_params'](settings)


def filename_for_scan(settings, year, descr):
    return "{}-{}-{}.hdf5".format(settings['id'], year, descr)


def filename_for_scan_summary(settings):
    return "{}-scan.hdf5".format(settings['id'])


def filename_for_forecast(settings, date, descr):
    date_str = date.strftime("%Y-%m-%d")
    return "{}-{}-{}.hdf5".format(settings['id'], date_str, descr)


def filename_for_cache(settings, year, descr):
    return "cache-{}-{}-{}.hdf5".format(settings['id'], year, descr)


def filename_for_forecast_json(settings, year, descr):
    json_file = "{}-{}-{}.json".format(settings['id'], year, descr)
    return os.path.join(settings['json_dir'], json_file)


def filename_for_scan_json(settings, year):
    json_file = "{}-{}.json".format(settings['id'], year)
    return os.path.join(settings['json_dir'], json_file)


def find_scan_files(settings, out_dir):
    regex = '^' + settings['id'] + r'-\d{4}-.*\.hdf5$'
    return [os.path.join(out_dir, f) for f in os.listdir(out_dir)
            if re.match(regex, f)]


def find_forecast_files(settings, descr, out_dir):
    # The location ID may be a prefix of another location, and so we need to
    # ensure that the text after the location ID is the forecasting date.
    regex = '^' + settings['id'] + r'-\d{4}-\d{2}-\d{2}-' + descr + '\\.hdf5$'
    return [os.path.join(out_dir, f) for f in os.listdir(out_dir)
            if re.match(regex, f)]


def override(name, default, cmdline_args, locn_settings):
    if cmdline_args[name] is not None:
        return cmdline_args[name]
    elif name in locn_settings:
        return locn_settings[name]
    else:
        return default


def om_values(values, year):
    """
    Return a list of values for an observation model parameter.

    :param values: Values may be scalar, a list of scalars, or a dictionary
        that maps years to a scalar or list of scalars.
    :param year: The calendar year of the simulation.

    :raises ValueError: If the values are a dictionary and the year is not a
        valid key, or if there are duplicate values.
    :raises TypeError: If the values are not one of the types defined above.
    """
    if isinstance(values, dict):
        # The values are stored in a dictionary and are indexed by year.
        if year in values:
            values = values[year]
        else:
            raise ValueError('Invalid year {}'.format(year))
    try:
        # Try treating the values as a list of values.
        values = list(values)
    except TypeError:
        # If that fails, assume we have a scalar value.
        return [values]
    # Detect duplicate values.
    if len(set(values)) == len(values):
        return values
    else:
        raise ValueError("Duplicate values in {}".format(values))


def om_param_names(settings, param_key):
    om_pnames = sorted(settings[param_key].keys())
    if len(om_pnames) == 0:
        raise ValueError("No parameters defined for '{}'".format(param_key))
    else:
        return om_pnames


def descr_for(settings, param_key, om_param_vals):
    om_pnames = om_param_names(settings, param_key)
    disp = [settings['om_name'][pn] for pn in om_pnames]
    fmt = [settings['om_format'][pn] for pn in om_pnames]

    out_fields = ['{{disp[{0}]}}-{{val[{0}]:{{fmt[{0}]}}}}'.format(i)
                  for i in range(len(om_pnames))]
    out_fmt = "-".join(out_fields)
    return out_fmt.format(disp=disp, val=om_param_vals, fmt=fmt)


def param_generator(settings, year, param_key, extra_arg=None):
    """
    Return a generator that performs a parameter sweep for observation models.
    This will yield tuples of the form ``(om_params, descr, extra)``, where
    the elements comprise:

    - ``om_params``: a dictionary that maps observation model parameter names
        to their values, which can be applied by:
        ``observer.define_params(params, **om_params)``.
    - ``descr``: a string that includes the name and value of each observation
        model parameter, for inclusion in the output file name.
    - ``extra``: the value of the ``extra_arg`` parameter.

    :param settings: The location-specific forecasting parameters.
    :param year: The calendar year (default is the current year).
    :param param_key: The key in ``settings`` that defines the sweep.
    :param extra_arg: Additional content to include in each tuple.
    """
    om_pnames = om_param_names(settings, param_key)
    disp = [settings['om_name'][pn] for pn in om_pnames]
    fmt = [settings['om_format'][pn] for pn in om_pnames]

    out_fields = ['{{disp[{0}]}}-{{val[{0}]:{{fmt[{0}]}}}}'.format(i)
                  for i in range(len(om_pnames))]
    out_fmt = "-".join(out_fields)

    om_scan = [om_values(settings[param_key][pn], year) for pn in om_pnames]
    for om_pvals in itertools.product(*om_scan):
        om_params = dict(zip(om_pnames, om_pvals))
        out_desc = out_fmt.format(disp=disp, val=om_pvals, fmt=fmt)
        yield (om_params, out_desc, extra_arg)


def locn_forecasts(settings, year=None, extra_arg=None):
    """
    Return a generator that performs a parameter sweep for live forecasts.
    This will yield tuples of the form ``(om_params, descr, extra)``, where
    the elements comprise:

    - ``om_params``: a dictionary that maps observation model parameter names
        to their values, which can be applied by:
        ``observer.define_params(params, **om_params)``.
    - ``descr``: a string that includes the name and value of each observation
        model parameter, for inclusion in the output file name.
    - ``extra``: the value of the ``extra_arg`` parameter.

    :param settings: The location-specific forecasting parameters.
    :param year: The calendar year (default is the current year).
    :param extra_arg: Additional content to include in each tuple.
    """
    if year is None:
        year = datetime.date.today().year

    return param_generator(settings, year, 'forecast', extra_arg)


def locn_scan(settings, year, extra_arg=None):
    """
    Return a generator that performs a parameter sweep for a single year in a
    single location.
    This will yield tuples of the form ``(om_params, descr, extra)``, where
    the elements comprise:

    - ``om_params``: a dictionary that maps observation model parameter names
        to their values, which can be applied by:
        ``observer.define_params(params, **om_params)``.
    - ``descr``: a string that includes the name and value of each observation
        model parameter, for inclusion in the output file name.
    - ``extra``: the value of the ``extra_arg`` parameter.

    :param settings: The location-specific forecasting parameters.
    :param year: The calendar year.
    :param extra_arg: Additional content to include in each tuple.

    :raises ValueError: If the year is not in ``settings['scan_years']``.
    """
    if year not in settings['scan_years']:
        raise ValueError("Invalid year {} for {}".format(
            year, settings['name']))

    return param_generator(settings, year, 'scan', extra_arg)


def multi_scan(locn_ids=None, years=None, extra_arg=None):
    """
    Return a generator that performs a parameter sweep over multiple locations
    and multiple years.
    It yields tuples of the form ``(om_params, descr, locn, year, extra)``,
    where the elements comprise:

    - ``om_params``: a dictionary that maps observation model parameter names
        to their values, which can be applied by:
        ``observer.define_params(params, **om_params)``.
    - ``descr``: a string that includes the name and value of each observation
        model parameter, for inclusion in the output file name.
    - ``locn``: the location ID.
    - ``year``: the calendar year.
    - ``extra``: the value of the ``extra_arg`` parameter.

    :param locn_ids: The location IDs over which to scan.
    :param years: The calendar years over which to scan.
    :param extra_arg: Additional content to include in each tuple.

    :raises ValueError: If the year is not in ``settings['scan_years']``.
    """
    valid_locns = locations()

    if locn_ids is None:
        locn_ids = sorted(valid_locns)

    invalid = [locn for locn in locn_ids if locn not in valid_locns]
    if invalid:
        raise ValueError("invalid locations: {}".format(", ".join(invalid)))

    for locn_id in locn_ids:
        s = local(locn_id)

        if years is None:
            l_years = s['scan_years']
        else:
            l_years = years

        invalid = [y for y in l_years if y not in s['scan_years']]
        if invalid:
            raise ValueError("invalid years {} for {}".format(years, locn_id))

        for year in l_years:
            gen = param_generator(s, year, 'scan', extra_arg)
            for (om_params, descr, extra) in gen:
                yield (om_params, descr, locn_id, year, extra)
