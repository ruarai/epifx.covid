#!/usr/bin/env python
"""
This script performs the following series of actions:

1. Downloads the latest Reff(t) trajectories and COVID-19 case data, as
   provided by Nick Golding via Dropbox shared folders.

2. Generates the input Reff(t) and COVID-19 case data files for each
   jurisdiction.

3. Generates a forecast for each jurisdiction.

4. Collects the forecast results and saves them to several CSV files:

   (a) Forecast credible intervals;

   (b) Parameter credible intervals;

   (c) Daily COVID-19 case counts;

   (d) Forecast sample trajectories; and

   (e) Parameter and state variable sample trajectories.

5. Generates the output file for the ensemble forecast:

      moss_backcast_and_forecast_samples_YYYY-MM-DD.csv

6. Generates validation plots.
"""

import argparse
import datetime
import epifx
import epifx.cmd.decl_fs
import h5py
import multiprocessing
import numpy as np
import os
import os.path
import re
import subprocess
import sys


def main(args=None):
    """
    Perform all of the actions described at the top of this file.
    """
    p = get_parser()
    options = p.parse_args(args)
    
    # Only generate live forecasts for the baseline model.
    toml_files = [options.ff]
    
    forecast_args = [arg for f in toml_files for arg in ['-c', f]]
    num_cpus = multiprocessing.cpu_count()
    if num_cpus >= 8:
        # Run each forecast in parallel if there are sufficient CPUs.
        forecast_args.extend(['--spawn', '8'])
        
    epifx.cmd.decl_fs.main(args=forecast_args)

    return 0


def get_parser():
    p = argparse.ArgumentParser()
        
    p.add_argument('--ff',
        action = 'store', type = str, help = 'The forecast file to read (default forecast.toml).',
        default = 'forecast.toml')

    return p


if __name__ == "__main__":
    sys.exit(main())
