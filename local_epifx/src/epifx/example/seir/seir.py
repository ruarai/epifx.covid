"""Example SEIR forecasts."""

from epifx.cmd.declarative import config_from_string
from epifx.cmd.decl_fs import forecast_iter, run

import datetime
import logging
import pkgutil


logging.basicConfig(level=logging.INFO)

toml_file = 'seir.toml'
toml_data = pkgutil.get_data('epifx.example.seir', toml_file).decode()
with open(toml_file, mode='w') as f:
    f.write(toml_data)

config = config_from_string(toml_data)

obs_file = 'weekly-cases.ssv'
obs_data = pkgutil.get_data('epifx.example.seir', obs_file).decode()
with open(obs_file, mode='w') as f:
    f.write(obs_data)

location_names = ['test']
forecast_from = datetime.datetime(2014, 4, 1)

forecasts = forecast_iter(config, location_names, forecast_from)
for forecast in forecasts:
    run(forecast)
