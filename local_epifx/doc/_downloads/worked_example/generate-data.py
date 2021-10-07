#!/usr/bin/env python

import epifxlocns

import datetime
import epifx.obs
import epifx.summary
import numpy as np
import pypfilt
import pypfilt.summary
import sys


def random_obs(rnd, om_params, exp_value):
    bg_var = om_params['bg_var']
    disp = om_params['disp']
    # Ensure that the variance is not smaller than the background variance.
    nb_var = exp_value + exp_value * exp_value / disp
    if nb_var < bg_var:
        disp = exp_value * exp_value / (bg_var - exp_value)
    nb_pr = disp / (disp + exp_value)
    return rnd.negative_binomial(disp, nb_pr)


def main(args=None):
    location = 'some-city'
    year = 2017
    settings = epifxlocns.local_settings(location)
    params = epifxlocns.get_locn_params(settings)
    params['hist']['px_count'] = 1
    params['hist']['wind_size'] = 0
    params['hist']['wind_shift'] = 0

    # Fix the model priors so that the precise ground truth is known.
    # Set R0 to 1.4.
    params['param_min'][4] = 1.4
    params['param_max'][4] = 1.4
    # Set the (inverse of) the incubation period.
    params['param_min'][5] = 1.0
    params['param_max'][5] = 1.0
    # Set the (inverse of) the infectious period.
    params['param_min'][6] = 0.5
    params['param_max'][6] = 0.5
    # Keep eta fixed at 1 (i.e., enforce homogeneous mixing).
    params['param_min'][7] = 1.0
    params['param_max'][7] = 1.0
    # No seasonal forcing.
    params['param_min'][8] = 0.0
    params['param_max'][8] = 0.0
    # Introduce the first infection after 15 weeks.
    params['param_min'][9] = 105
    params['param_max'][9] = 105
    # Update the priors.
    params['prior'] = params['model'].priors(params)

    start = epifxlocns.get_start_date(year)
    until = epifxlocns.get_until_date(year)
    fake_obs = {
        'date': datetime.datetime(2017, 1, 1),
        'unit': 'Weekly Cases',
        'period': 7,
        'value': 0,
    }
    obs_model = epifx.obs.PopnCounts('Weekly Cases', 7)
    bg_obs = 5
    bg_var = 5
    pr_obs = 0.01
    disp = 100
    obs_model.define_params(params, bg_obs, pr_obs, disp, bg_var)
    om_params = params['obs']['Weekly Cases']
    summary = pypfilt.summary.HDF5(params, obs_list=[fake_obs])
    peak_monitor = epifx.summary.PeakMonitor()
    summary.add_tables(epifx.summary.ExpectedObs(peak_monitor, probs=[0]))
    state = pypfilt.run(params, start, until, [], summary)
    exp_obs = state['summary']['expected_obs']
    rnd = np.random.RandomState(year)
    with open('case-counts.ssv', 'w') as f:
        f.write("year date       cases\n")
        for row in exp_obs[::7][1:]:
            date_bs = row['date']
            date = params['time'].from_dtype(date_bs).date()
            exp_cases = row['ymin']
            rnd_cases = random_obs(rnd, om_params, exp_cases)
            f.write("{} {} {:.0f}\n".format(year, date, rnd_cases))
    return 0


if __name__ == "__main__":
    sys.exit(main())
