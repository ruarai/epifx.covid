"""
Provides a declarative means of defining forecasts.

The purpose of this module is to allow users to define and run forecasting
simulations **without writing any Python code**, by instead defining all of
the necessary settings and parameters in a `TOML`_ file.
"""

import copy
import logging
import importlib
import os
import sys
from pathlib import Path
from typing import Any, Callable, Dict, NamedTuple, Optional
import toml

from .context import Context
from .params import default_params
from . import io
from . import obs
from . import summary


ParamsFn = Callable[[], dict]
"""Functions that return simulation parameters."""


SummaryFn = Callable[..., summary.HDF5]
"""Functions that define summary statistics."""


class ObsModel(NamedTuple):
    """
    Observation model settings.

    :param obs.Obs obs_model: The observation model.
    :param str obs_file: The name of the file that contains the observations.
    :param file_args: Keyword arguments for the observation model's
        ``from_file()`` method.
    :type file_args: Dict[str, Any]
    :param parameters: The observation model parameter values.
    :type parameters: Dict[str, float]
    :param names: The name to use for each parameter in output file names.
    :type names: Dict[str, str]
    :param formats: The format specifier to use for each parameter value in
        output file names.
    :type formats: Dict[str, str]
    :param dict plot_settings: Settings for interactive forecast plots.
    """
    obs_model: obs.Obs
    obs_file: str
    file_args: Dict[str, Any]
    parameters: Dict[str, float]
    names: Dict[str, str]
    formats: Dict[str, str]
    plot_settings: Dict


class Scenario(NamedTuple):
    """
    Scenario-specific settings.

    :param str name: The scenario name.
    :param params_fn: A function that returns the simulation parameters.
    :type params_fn: :py:data:`ParamsFn`
    :param summary_fn: A function that defines the summary statistics.
    :type summary_fn: :py:data:`SummaryFn`
    :param summary_args: Arguments for the summary statistic function.
    :type summary_args: Dict[str, Any]
    :param observation_models: The observation model settings.
    :type observation_models: Dict[str, pypfilt.config.ObsModel]
    """
    name: str
    params_fn: ParamsFn
    summary_fn: SummaryFn
    summary_args: Dict[str, Any]
    observation_models: Dict[str, ObsModel]


class Config(NamedTuple):
    """
    A configuration.

    :param scenarios: Scenario-specific settings.
    :type scenarios: Dict[str, Scenario]
    :param str config_string: The configuration file contents.
    :param dict config_toml: The parsed configuration file contents.
    """
    scenarios: Dict[str, Scenario]
    config_string: str
    config_toml: Dict


class _Settings(NamedTuple):
    """
    The configuration settings that apply in a particular context, such as
    global settings or scenario-specific settings.
    """
    time: Any
    model: Any
    summary_fn: Any
    params: Dict[str, Any]
    priors: Dict[str, Any]
    bounds: Dict[str, Any]
    summary_args: Dict[str, Any]
    summary_monitors: Dict[str, Any]
    summary_tables: Dict[str, Any]


def from_file(toml_path: str) -> Config:
    """Read a forecast configuration from a TOML file."""
    with open(toml_path, encoding='utf-8') as f:
        toml_string = f.read()
    return from_string(toml_string)


def from_string(toml_string: str) -> Config:
    """Read a forecast configuration from a string."""
    cfg_data = toml.loads(toml_string)

    scenarios_data = __get(cfg_data, 'scenario')
    defaults = {k: v for (k, v) in cfg_data.items() if k != 'scenario'}

    scenarios = {
        scen_id: make_scenario(scen_id, scen_data, copy.deepcopy(defaults))
        for (scen_id, scen_data) in scenarios_data.items()
    }

    config = Config(
        scenarios=scenarios,
        config_string=toml_string,
        config_toml=cfg_data,
    )
    return config


def make_settings(cfg_data: Dict[str, Any]) -> _Settings:
    """Return the settings for a given configuration."""
    logger = logging.getLogger(__name__)

    if 'component' in cfg_data and 'components' not in cfg_data:
        logger.warning('"%s": did you mean "%s"?',
                       'component', 'components')

    components = __get(cfg_data, 'components')
    model = instantiate(__get(components, 'model', 'Components'))
    time = instantiate(__get(components, 'time', 'Components'))
    summary_fn = lookup(__get(components, 'summary', 'Components'))
    for key in components:
        if key not in ['model', 'time', 'summary']:
            logger.warning('unrecognised component "%s"', key)

    model_priors = {}
    if 'model' in cfg_data:
        if 'priors' in cfg_data['model']:
            model_priors = cfg_data['model']['priors']
        elif 'prior' in cfg_data['model']:
            logger.warning('"%s": did you mean "%s"?',
                           'model.prior', 'model.prior')

    model_bounds = {}
    if 'model' in cfg_data:
        if 'bounds' in cfg_data['model']:
            model_bounds = cfg_data['model']['bounds']
        elif 'bound' in cfg_data['model']:
            logger.warning('"%s": did you mean "%s"?',
                           'model.bound', 'model.bounds')

    summary_args = {}
    monitors = {}
    tables = {}
    if 'summary' in cfg_data:
        if 'init' in cfg_data['summary']:
            summary_args = cfg_data['summary']['init']
        if 'monitors' in cfg_data['summary']:
            for (name, data) in cfg_data['summary']['monitors'].items():
                mon_name = __get(data, 'model',
                                 'summary.monitors.{}'.format(name))
                mon_args = data.get('init', {})
                monitors[name] = instantiate(mon_name, **mon_args)
        elif 'monitor' in cfg_data['summary']:
            logger.warning('"%s": did you mean "%s"?',
                           'summary.monitor', 'summary.monitors')
        if 'tables' in cfg_data['summary']:
            for (name, data) in cfg_data['summary']['tables'].items():
                tbl_name = __get(data, 'model',
                                 'summary.tables.{}'.format(name))
                tbl_args = data.get('init', {})
                tables[name] = instantiate(tbl_name, **tbl_args)
        elif 'table' in cfg_data['summary']:
            logger.warning('"%s": did you mean "%s"?',
                           'summary.table', 'summary.tables')

    param_data = __get(cfg_data, 'parameters')

    return _Settings(
        time=time,
        model=model,
        summary_fn=summary_fn,
        params=param_data,
        priors=model_priors,
        bounds=model_bounds,
        summary_args=summary_args,
        summary_monitors=monitors,
        summary_tables=tables)


def __get(data: Dict[str, Any], key: str, src: Optional[str] = None) -> Any:
    """
    Get a value from a dictionary; if the key does not exist, raise an
    exception that identifies the missing key and the configuration section in
    which it was expected.

    :param data: The dictionary from which to get the value.
    :param key: The key name.
    :param src: The configuration section associated with the dictionary.
    """
    try:
        return data[key]
    except KeyError:
        if src is None:
            src = 'Configuration'
        if src:
            raise ValueError('{}: "{}" is missing'.format(src, key))
        else:
            raise ValueError('"{}" is missing'.format(key))


def override_dict(defaults, overrides):
    """
    Override a dictionary with values in another dictionary. This will
    recursively descend into matching nested dictionaries.

    :param dict defaults: The original values.
    :param dict overrides: The overriding values.
    :return: The modified ``params`` dictionary.
    :rtype: Dict
    """
    for (key, value) in overrides.items():
        if isinstance(value, dict):
            # NOTE: if the parameter already exists, we can only use it as a
            # default value if it is a dictionary.
            if key in defaults and isinstance(defaults[key], dict):
                sub_defaults = defaults[key]
            else:
                sub_defaults = {}
            defaults[key] = override_dict(sub_defaults, value)
        else:
            defaults[key] = value
    return defaults


def validate_priors(settings: _Settings) -> None:
    """
    Ensure each prior comprises a function name and arguments dictionary,
    otherwise raise a ValueError.

    Note that this doesn't enforce that each prior corresponds to a known
    model parameter, because this would prevent us from supporting prior
    distributions that are expressed in terms of **transformed** parameters
    (such as reciprocals of rate parameters).
    """
    logger = logging.getLogger(__name__)
    for (name, info) in settings.priors.items():
        if 'function' not in info:
            raise ValueError('Missing prior function for {}'.format(name))
        elif not isinstance(info['function'], str):
            raise ValueError('Invalid prior function for {}'.format(name))
        if 'args' not in info:
            raise ValueError('Missing prior arguments for {}'.format(name))
        elif not isinstance(info['args'], dict):
            raise ValueError('Invalid prior arguments for {}'.format(name))
        if len(info) != 2:
            extra_keys = [k for k in info if k not in ['function', 'args']]
            logger.warning('Extra prior keys for %s: %s', name, extra_keys)


def validate_bounds(settings: _Settings) -> None:
    """Ensure bounds only refer to known parameters, or raise a ValueError."""
    descr = {name: ix
             for (ix, (name, _smooth, _lower, _upper)) in
             enumerate(settings.model.describe())}
    unknown = set(settings.bounds.keys()) - set(descr.keys())
    if unknown:
        raise ValueError('Unknown model parameters: {}'.format(unknown))
    for (name, info) in settings.bounds.items():
        if 'min' not in info or 'max' not in info:
            raise ValueError('No parameter bounds for model parameter {}'
                             .format(name))


def make_obs_models(scen_data: dict, scen_name: str) -> Dict[str, ObsModel]:
    """Create the observation models for a scenario."""
    if 'observations' not in scen_data:
        raise ValueError('No observation models for {}'.format(scen_name))

    obs_models = {}
    for (obs_unit, obs_data) in scen_data['observations'].items():
        source = 'observations.{}'.format(obs_unit)
        kwargs = {'obs_unit': obs_unit}
        if 'init' in obs_data:
            for (name, value) in obs_data['init'].items():
                kwargs[name] = value

        om = instantiate(__get(obs_data, 'model', source), **kwargs)

        plot_settings = {
            'scenario_label': scen_name,
            'axis_label': '',
            'axis_precision': 0,
            'datum_label': '',
            'datum_precision': 0,
        }
        for (name, value) in obs_data.get('plot', {}).items():
            plot_settings[name] = value

        obs_models[obs_unit] = ObsModel(
            obs_model=om,
            obs_file=__get(obs_data, 'file', source),
            file_args=obs_data.get('file_args', {}),
            parameters=__get(obs_data, 'parameters', source),
            formats=__get(obs_data, 'format', source),
            names=__get(obs_data, 'name', source),
            plot_settings=plot_settings,
        )

    return obs_models


def make_column_init_fn(num_values: int) -> Callable[[Context, int], None]:
    """
    Return a function that will associate each particle with a single column
    in the named lookup table.
    """

    def init_fn(ctx, col_vec):
        rnd = ctx.component['random']['hist_extra_cols']
        col_vec[:] = rnd.integers(num_values, size=col_vec.shape)

    return init_fn


def create_lookup_tables(scen_data, settings: _Settings, params) -> None:
    """Create all of the lookup tables and add them to the parameters."""
    if 'lookup_tables' in scen_data:
        data_dir = Path(__get(params, 'data_dir', 'Default parameters'))
        for (name, filename) in scen_data['lookup_tables'].items():
            data = io.read_lookup_table(data_dir / filename,
                                        settings.time)
            table = io.Lookup(data)
            params['data']['lookup'][name] = data
            params['component']['lookup'][name] = table

            # Check whether a column should be added to the state matrix
            # so that particles can draw samples from this lookup table.
            # NOTE: this feature is used by epifx to allow observation model
            # parameters to be sampled independently for each particle.
            if 'sample_lookup_tables' in scen_data:
                if name in scen_data['sample_lookup_tables']:
                    num_vals = table.value_count()
                    init_fn = make_column_init_fn(num_vals)
                    params['hist']['extra_col_fns'][name] = init_fn
                    params['hist']['extra_cols'] += 1
                    if 'sample_lookup_tables' not in params:
                        params['sample_lookup_tables'] = {}
                    params['sample_lookup_tables'][name] = (
                        params['hist']['extra_cols'] - 3)

    # Ensure all tables associated with a state column exist.
    if 'sample_lookup_tables' in scen_data:
        for name in scen_data['sample_lookup_tables']:
            if name not in params['component']['lookup']:
                raise ValueError('Missing table: {}'.format(name))


def define_summary_tables(settings, params):
    """Add summary tables and monitors to the parameters."""
    for (name, monitor) in settings.summary_monitors.items():
        params['component']['summary_monitor'][name] = monitor

    for (name, table) in settings.summary_tables.items():
        params['component']['summary_table'][name] = table


def make_prior_fn(fn_name, args):
    """
    Return a function that draws samples from the specified distribution.

    :param fn_name: The name of a ``numpy.random.Generator`` method used to
        generate samples.
    :param args: A dictionary of keyword arguments.

    As a special case, ``fn_name`` may be set to ``'inverse_uniform'`` to
    sample from a uniform distribution and then take the reciprocal:

    .. math:: X \\sim \\frac{1}{\\mathcal{U}(a, b)}

    The bounds ``a`` and ``b`` may be specified by the following keyword
    arguments:

    + :code:`a = args['low']` **or** :code:`a = 1 / args['inv_low']`
    + :code:`b = args['high']` **or** :code:`b = 1 / args['inv_high']`
    """
    if fn_name == 'inverse_uniform':
        if 'low' in args:
            low = args['low']
        else:
            low = 1 / args['inv_low']

        if 'high' in args:
            high = args['high']
        else:
            high = 1 / args['inv_high']

        return lambda r, size=None: 1 / r.uniform(
            low=low, high=high, size=size)
    else:
        return lambda r, size=None: getattr(r, fn_name)(**args, size=size)


def make_params_fn(scen_data, settings: _Settings,
                   scen_id, scen_name) -> Callable[[], dict]:
    """
    Return a no-argument function that will return the simulation parameters.
    """
    px_count = __get(settings.params, 'particles', 'Parameters')
    max_days = __get(settings.params, 'max_days', 'Parameters')
    prng_seed = __get(settings.params, 'prng_seed', 'Parameters')
    ignore = ['particles', 'prng_seed', 'max_days']
    overrides = {k: v for (k, v) in settings.params.items()
                 if k not in ignore}

    descr = {name: ix
             for (ix, (name, _smooth, _lower, _upper)) in
             enumerate(settings.model.describe())}

    def params_fn():
        params = default_params(settings.model,
                                settings.time,
                                max_days=max_days,
                                px_count=px_count,
                                prng_seed=prng_seed)
        params = override_dict(params, overrides)

        # NOTE: ensure the start and end of the simulation period are
        # correctly defined. If they are strings, use time.from_unicode() to
        # convert them to time values, otherwise ensure that they are valid
        # time values.
        if isinstance(params['time']['start'], str):
            params['time']['start'] = settings.time.from_unicode(
                params['time']['start'])
        elif not settings.time.is_instance(params['time']['start']):
            raise ValueError('Invalid value for time.start: "{}"'.format(
                type(params['time']['start'])))

        if isinstance(params['time']['until'], str):
            params['time']['until'] = settings.time.from_unicode(
                params['time']['until'])
        elif not settings.time.is_instance(params['time']['until']):
            raise ValueError('Invalid value for time.until: "{}"'.format(
                type(params['time']['until'])))

        if settings.bounds:
            for (name, info) in settings.bounds.items():
                # NOTE: apply parameter bounds.
                param_ix = descr[name]
                params['model']['param_min'][param_ix] = info['min']
                params['model']['param_max'][param_ix] = info['max']

        # Allow priors to specify an RNG sampling function and arguments.
        if settings.priors:
            for (name, prior) in settings.priors.items():
                params['model']['prior'][name] = make_prior_fn(
                    prior['function'], prior['args'])

        create_lookup_tables(scen_data, settings, params)
        define_summary_tables(settings, params)

        # Record the scenario ID and name in the parameters dictionary, so
        # that it is recorded in each output file.
        params['scenario'] = {
            'id': scen_id,
            'name': scen_name,
        }

        return params

    return params_fn


def make_scenario(scen_id: str, scen_data: dict, defaults: dict) -> Scenario:
    """
    Return the settings for a single scenario.
    """
    scen_name = __get(scen_data, 'name', 'scenario: {}'.format(scen_id))
    cfg_data = override_dict(defaults, scen_data)
    settings = make_settings(cfg_data)
    obs_models = make_obs_models(cfg_data, scen_name)
    validate_priors(settings)
    validate_bounds(settings)
    params_fn = make_params_fn(cfg_data, settings, scen_id, scen_name)

    return Scenario(
        name=scen_name,
        params_fn=params_fn,
        summary_fn=settings.summary_fn,
        summary_args=settings.summary_args,
        observation_models=obs_models)


def lookup(full_name):
    """
    Retrieve an object from a fully-qualified name.

    :param str full_name: The fully-qualified name.

    :Examples:

    >>> summary_fn = lookup('pypfilt.summary.HDF5')
    """
    last_dot = full_name.rfind('.')
    if last_dot < 0:
        raise ValueError('No module name in "{}"'.format(full_name))
    module_name = full_name[:last_dot]
    value_name = full_name[last_dot + 1:]

    try:
        module = importlib.import_module(module_name)
    except ModuleNotFoundError:
        logger = logging.getLogger(__name__)
        logger.debug('Could not import "{}"'.format(module_name))
        logger.debug('Trying again in the working directory ...')
        cwd = os.getcwd()
        sys.path.append(cwd)
        try:
            module = importlib.import_module(module_name)
        finally:
            sys.path.pop()
        logger.debug('Successfully imported "{}"'.format(module_name))

    # NOTE: will raise AttributeError if the attribute does not exist.
    value = getattr(module, value_name)
    return value


def instantiate(full_name, *args, **kwargs):
    """
    Instantiate an object from a class name.

    :param str full_name: The fully-qualified class name.
    :param \\*args: Positional constructor arguments (optional).
    :param \\**kwargs: Named constructor arguments (optional).

    :Examples:

    >>> time = instantiate('pypfilt.Datetime', date_fmt='%Y-%m-%d')
    """
    object_class = lookup(full_name)
    if not callable(object_class):
        raise ValueError('The value "{}" is not callable'.format(full_name))
    try:
        return object_class(*args, **kwargs)
    except TypeError:
        logger = logging.getLogger(__name__)
        logger.error('Attempted to call "{}" with arguments:'
                     .format(full_name))
        for arg in args:
            print('    {}'.format(arg))
        for (name, value) in kwargs.items():
            print('    {} = {}'.format(name, value))
        raise
