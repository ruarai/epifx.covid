"""Sweep over forecast settings."""

import itertools
from pathlib import Path
from typing import Dict, List, NamedTuple, Sequence

from . import config


class Forecasts(NamedTuple):
    """
    Forecast settings.

    :param str scenario_id: The scenario identifier.
    :param str scenario_name: A descriptive name for the scenario.
    :param dict params: The simulation parameters.
    :param List[dict] all_observations: All of the available observations.
    :param List[Sequence[dict]]: observation_streams: The observations,
        separated into streams for each observation model.
    :param str file_descriptor: A descriptor to identify files associated with
        these forecasts.
    :param List[str] config_strings: The contents of the configuration files
        in the forecasting sweep.
    """
    scenario_id: str
    scenario_name: str
    params: dict
    observation_models: Dict[str, config.ObsModel]
    all_observations: List[dict]
    observation_streams: List[Sequence[dict]]
    file_descriptor: str
    config_strings: List[str]


def forecasts(configs, load_obs=True):
    """
    Iterate over the forecast settings for each scenario.

    :param configs: The forecasting configuration(s).
    :type configs: Union[pypfilt.config.Config, Sequence[pypfilt.config.Config]]
    :param bool load_obs: Whether to load observations from disk.
    :rtype: Iterator[Forecasts]
    """
    if isinstance(configs, config.Config):
        configs = [configs]

    config_strings = [cfg.config_string for cfg in configs]

    for cfg in configs:
        for (scenario_id, scenario_config) in cfg.scenarios.items():
            for (obs_models, descr) in obs_model_combs(scenario_config):
                (params, all_obs, obs_streams) = get_params(scenario_config,
                                                            obs_models,
                                                            load_obs=load_obs)
                params['component']['summary'] = scenario_config.summary_fn(
                    params, all_obs, **scenario_config.summary_args)
                forecasts = Forecasts(
                    scenario_id=scenario_id,
                    scenario_name=scenario_config.name,
                    params=params,
                    observation_models=scenario_config.observation_models,
                    all_observations=all_obs,
                    observation_streams=obs_streams,
                    file_descriptor=descr,
                    config_strings=config_strings)
                yield forecasts


def forecasts_mp(configs):
    """
    Iterate over the forecast settings for each scenario, yielding tuples
    that can be serialised and sent to other Python processes.
    Each process can then retrieve the associated :class:`Forecasts` value by
    calling :func:`get_forecasts_mp`.

    :param configs: The forecasting configuration(s).
    :type configs: Union[pypfilt.config.Config, Sequence[pypfilt.config.Config]]
    """
    for (ix, forecast) in enumerate(forecasts(configs, load_obs=False)):
        yield (ix, forecast.config_strings, forecast.scenario_id,
               forecast.file_descriptor)


def get_forecasts_mp(mp_tuple, load_obs=True):
    """
    Return the forecast settings for a tuple yielded by :func:`forecasts_mp`.

    :param mp_tuple: A value returned by :func:`forecasts_mp`.
    :param bool load_obs: Whether to load observations from disk.
    :rtype: Forecasts
    """
    (ix, config_strings, scenario_id, file_descriptor) = mp_tuple
    configs = [config.from_string(s) for s in config_strings]
    for (i, forecast) in enumerate(forecasts(configs, load_obs=load_obs)):
        if i == ix:
            # NOTE: ensure that the details are consistent with the original
            # Forecasts object generated by the original forecasts iterator.
            if forecast.scenario_id != scenario_id:
                raise ValueError('Forecast #{} expected {} found {}'.format(
                    ix, scenario_id, forecast.scenario_id))
            if forecast.file_descriptor != file_descriptor:
                raise ValueError('Forecast #{} expected {} found {}'.format(
                    ix, file_descriptor, forecast.file_descriptor))
            return forecast
    raise ValueError('Could not find forecast #{}'.format(ix))


def obs_model_combs(config):
    """
    Iterate over every combination of observation model parameters for each of
    the defined observation models.

    :param config: The forecasting configuration.
    :type config: declarative.Config
    """
    obs_models = [
        obs_model_parameters(obs_unit, obs_config)
        for (obs_unit, obs_config) in config.observation_models.items()]
    for obs_models in itertools.product(*obs_models):
        # Join the descriptors for each observation model.
        descrs = [descr for (_obs_unit, _params, descr) in obs_models]
        descr = '-'.join(descrs)
        # Collect all of the observation model parameters.
        obs_config = {obs_unit: params
                      for (obs_unit, params, _) in obs_models}
        yield(obs_config, descr)


def obs_model_parameters(obs_unit, obs_config):
    """
    Iterate over every combination of observation model parameters for a
    single observation model.

    :param str obs_unit: The observation units.
    :param dict obs_config: The observation model parameter values.
    """
    param_values = obs_config.parameters
    param_formats = obs_config.formats
    # NOTE: sort parameter names to enforce a consistent ordering.
    names = sorted(param_values.keys())
    display_names = obs_config.names

    out_fields = []
    for (ix, name) in enumerate(names):
        # NOTE: produce format strings such as 'bg-{val[0]:{fmt[bg_obs]}}'.
        field = '{0}-{{values[{1}]:{{formats[{2}]}}}}'.format(
            display_names[name], ix, name)
        out_fields.append(field)
    out_fmt = '-'.join(out_fields)

    # NOTE: the parameters must be scanned in their listed order.
    param_scan = [obs_model_parameter_values(param_values[name])
                  for name in names]
    for parameter_values in itertools.product(*param_scan):
        params_dict = dict(zip(names, parameter_values))
        description = out_fmt.format(values=parameter_values,
                                     formats=param_formats)
        yield (obs_unit, params_dict, description)


def obs_model_parameter_values(values):
    """
    Return a list containing all of the values for a single parameter.

    :param values: The values for the parameter.
    :type values: Union[List[Any], Any]
    """
    if isinstance(values, list):
        return values
    else:
        return [values]


def get_params(scenario_config, obs_config, load_obs=True):
    """
    Return the simulation parameters for the given combination of observation
    model parameters.

    :param scenario_config: The scenario configuration settings.
    :type scenario_config: pypfilt.config.Scenario
    :param obs_config: The parameter values for each observation model.
    :type obs_config: Dict[Dict[str, Any]]
    :param bool load_obs: Whether to load observations from disk.
    """
    params = scenario_config.params_fn()
    obs_streams = []

    for (obs_unit, obs_params) in obs_config.items():
        # Retrieve the observation model.
        obs_settings = scenario_config.observation_models[obs_unit]
        obs_model = obs_settings.obs_model
        # Apply the observation model parameters.
        params['obs'][obs_unit] = {}
        for (name, value) in obs_params.items():
            params['obs'][obs_unit][name] = value
        params['component']['obs'][obs_unit] = obs_model
        # Read the observations for this observation model.
        if load_obs:
            obs_file = Path(params['data_dir']) / obs_settings.obs_file
            # NOTE: also pass the time scale, before additional arguments.
            obs, tbl = obs_model.from_file(obs_file,
                                           params['component']['time'],
                                           **obs_settings.file_args)
            obs_streams.append(obs)
            params['data']['obs'][obs_unit] = tbl

    all_obs = list(itertools.chain.from_iterable(obs_streams))

    return (params, all_obs, obs_streams)