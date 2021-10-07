"""Test cases for the SEIR forecasting example."""

import datetime
import epifx.cmd.decl_fs as fs
import logging
import numpy as np
import os
import pkgutil
import pypfilt.config
import pypfilt.sweep


def two_forecast_dates(all_obs, fs_from):
    """Select only two forecasting dates, to reduce computation time."""
    first_obs = min(obs['date'] for obs in all_obs
                    if obs['date'] >= fs_from)
    twelve_weeks_later = first_obs + datetime.timedelta(days=12 * 7)
    return [first_obs, twelve_weeks_later]


def two_forecast_times(all_obs, fs_from):
    """Select only two forecasting times, to reduce computation time."""
    first_obs = min(obs['date'] for obs in all_obs
                    if obs['date'] >= fs_from)
    twelve_weeks_later = first_obs + 12 * 7
    return [first_obs, twelve_weeks_later]


def simulate_seir_observations():
    """Generate synthetic observations from a known model."""
    toml_file = 'seir.toml'
    pr_file = 'pr-obs.ssv'
    obs_file_dt = 'simulated-weekly-cases-datetime.ssv'
    obs_file_sc = 'simulated-weekly-cases-scalar.ssv'

    toml_data = pkgutil.get_data('epifx.example.seir', toml_file).decode()
    config = pypfilt.config.from_string(toml_data)

    pr_data = pkgutil.get_data('epifx.example.seir', pr_file).decode()
    with open(pr_file, mode='w') as f:
        f.write(pr_data)

    forecasts = list(pypfilt.sweep.forecasts(config, load_obs=False))
    assert len(forecasts) == 1
    forecast = forecasts[0]

    # NOTE: define the fixed ground truth for the model simulation.
    params = forecast.params
    params['model']['prior'] = {
        'R0': lambda r, size=None: r.uniform(low=1.45, high=1.45,
                                             size=size),
        'sigma': lambda r, size=None: r.uniform(low=0.25, high=0.25,
                                                size=size),
        'gamma': lambda r, size=None: r.uniform(low=0.25, high=0.25,
                                                size=size),
        'eta': lambda r, size=None: r.uniform(low=1.0, high=1.0,
                                                size=size),
        'alpha': lambda r, size=None: r.uniform(low=0.0, high=0.0,
                                                size=size),
        't0': lambda r, size=None: r.uniform(low=14.0, high=14.0,
                                             size=size),
    }
    obs_table = pypfilt.simulate_from_model(params, px_count=1)

    # Extract weekly observations from the simulated data.
    def to_date(bs):
        return datetime.datetime.strptime(bs.decode(),
                                          '%Y-%m-%d %H:%M:%S').date()

    obs_list = [(to_date(row['date']), row['value'].astype(int))
                for row in obs_table
                if to_date(row['date']).isoweekday() == 7]

    # Save date-indexed observations to disk.
    dt_obs = [(row[0].strftime('%Y-%m-%d'), row[1])
              for row in obs_list]
    dt_obs = np.array(dt_obs, dtype=[('date', 'O'), ('value', np.int_)])
    np.savetxt(obs_file_dt, dt_obs, fmt='%s %d',
               header='date cases', comments='')

    # Save day-indexed observations to disk.
    sc_obs = [(int(row[0].strftime('%-j')), row[1])
              for row in obs_list]
    sc_obs = np.array(sc_obs, dtype=[('day', np.int_), ('value', np.int_)])
    np.savetxt(obs_file_sc, sc_obs, fmt='%d %d',
               header='day cases', comments='')

    return (obs_list, obs_file_dt, obs_file_sc)


def test_simulate():
    """
    Generate synthetic observations from a known model, and check that the
    serialised results are consistent.
    """
    (obs_list, obs_file_dt, obs_file_sc) = simulate_seir_observations()

    peak_size = max(o[1] for o in obs_list)
    peak_time = [o[0] for o in obs_list if o[1] == peak_size][0]
    assert peak_size == 2678
    assert peak_time == datetime.date(2014, 9, 14)
    peak_day = int(peak_time.strftime('%-j'))

    # Check that the date-indexed peak is consistent with the above results.
    dt_cols = [pypfilt.io.date_column('date'), ('cases', int)]
    dt_obs = pypfilt.io.read_table(obs_file_dt, dt_cols)
    dt_mask = dt_obs['cases'] == peak_size
    assert np.sum(dt_mask) == 1
    assert dt_obs['date'][dt_mask].item().date() == peak_time

    # Check that the day-indexed peak is consistent with the above results.
    sc_cols = [('day', int), ('cases', int)]
    sc_obs = pypfilt.io.read_table(obs_file_sc, sc_cols)
    sc_mask = sc_obs['cases'] == peak_size
    assert np.sum(sc_mask) == 1
    assert sc_obs['day'][sc_mask] == peak_day

    # Check that the observations are the same.
    assert np.array_equal(sc_obs['cases'], dt_obs['cases'])

    # Clean up: remove created files.
    os.remove(obs_file_dt)
    os.remove(obs_file_sc)


def test_seeiir_forecast():
    """
    Use the SEEIIR forecasting example to compare peak size and time
    predictions at two forecasting dates.

    Note that the observation probability is set to 0.5 (much too high) and so
    we should only obtain sensible forecasts if the observation model is able
    to use the lookup table and obtain observation probabilities from the
    ``pr-obs.ssv`` data file.
    """
    logging.basicConfig(level=logging.INFO)

    toml_file = 'seeiir.toml'
    obs_file = 'weekly-cases.ssv'
    pr_file = 'pr-obs.ssv'

    toml_data = pkgutil.get_data('epifx.example.seir', toml_file).decode()
    config = pypfilt.config.from_string(toml_data)

    obs_data = pkgutil.get_data('epifx.example.seir', obs_file).decode()
    with open(obs_file, mode='w') as f:
        f.write(obs_data)

    pr_data = pkgutil.get_data('epifx.example.seir', pr_file).decode()
    with open(pr_file, mode='w') as f:
        f.write(pr_data)

    forecast_from = datetime.datetime(2014, 4, 1)

    # Check that there is only one set of forecasts (i.e., only one location
    # and only one set of observation model parameters).
    forecasts = pypfilt.sweep.forecasts(config)
    forecasts = list(forecasts)
    assert len(forecasts) == 1

    # Check that forecasts were run for two forecasting dates.
    forecast = forecasts[0]
    forecast_dates = two_forecast_dates(forecast.all_observations,
                                        forecast_from)
    state = fs.run(forecast, forecast_dates)
    fs_dates = list(state.keys())
    assert len(fs_dates) == 2

    fs_date_n1 = fs_dates[0]
    fs_date_n2 = fs_dates[1]

    # Retrieve the list of observations
    obs = state[fs_date_n1]['obs']
    peak_size = max(o['value'] for o in obs)
    peak_date = [o['date'] for o in obs if o['value'] == peak_size][0]

    # Check that the peak size and date is as expected.
    assert peak_size == 2678
    assert peak_date == datetime.datetime(2014, 9, 14)

    # Compare the forecast predictions to the observed peak size and date.
    forecast_n1 = state[fs_date_n1][fs_date_n1]['summary']
    forecast_n2 = state[fs_date_n2][fs_date_n2]['summary']
    dt_format = '%Y-%m-%d %H:%M:%S'

    # Ensure that all of the expected tables have been created, and that no
    # other tables have been created.
    expected_tables = {
        'model_cints', 'param_covar', 'pr_epi', 'forecasts', 'obs_llhd',
        'peak_size_acc', 'peak_time_acc', 'peak_cints', 'peak_ensemble',
        'obs/cases', 'exceed_500', 'exceed_1000', 'expected_obs'}
    tables_n1 = set(forecast_n1.keys())
    tables_n2 = set(forecast_n2.keys())
    assert tables_n1 == expected_tables
    assert tables_n2 == expected_tables
    # Ensure that no tables are empty.
    for name in expected_tables:
        shape_n1 = forecast_n1[name].shape
        shape_n2 = forecast_n1[name].shape
        assert len(shape_n1) == 1
        assert len(shape_n2) == 1
        assert shape_n1[0] > 0
        assert shape_n2[0] > 0
    # Ensure that the exceed_500 and exceed_1000 tables differ.
    pr_exc_low_n1 = forecast_n1['exceed_500'][()]['prob']
    pr_exc_high_n1 = forecast_n1['exceed_1000'][()]['prob']
    assert pr_exc_low_n1.shape == pr_exc_high_n1.shape
    assert not np.allclose(pr_exc_low_n1, pr_exc_high_n1)
    # Ensure that the cumulative probability of exceeding 500 cases is greater
    # than that of exceeding 1000 cases, until they both equal 1.0.
    cum_pr_low_n1 = np.cumsum(pr_exc_low_n1)
    cum_pr_high_n1 = np.cumsum(pr_exc_high_n1)
    mask_lt_1 = np.logical_and(cum_pr_low_n1 < 1.0, cum_pr_high_n1 < 1.0)
    mask_gt_0 = np.logical_or(cum_pr_low_n1 > 0.0, cum_pr_high_n1 > 0.0)
    mask = np.logical_and(mask_lt_1, mask_gt_0)
    assert np.all(cum_pr_low_n1[mask] > cum_pr_high_n1[mask])

    # The earlier forecast should include the peak size and time in its 95%
    # credible intervals.
    cints_n1 = forecast_n1['peak_cints']
    ci_n1 = cints_n1[cints_n1['prob'] == 95]
    ci_n1_size_lower = ci_n1['sizemin'].item()
    ci_n1_size_upper = ci_n1['sizemax'].item()
    ci_n1_date_lower = datetime.datetime.strptime(
        ci_n1['timemin'].item().decode(),
        dt_format)
    ci_n1_date_upper = datetime.datetime.strptime(
        ci_n1['timemax'].item().decode(),
        dt_format)
    assert ci_n1_size_lower <= peak_size <= ci_n1_size_upper
    assert ci_n1_date_lower <= peak_date <= ci_n1_date_upper

    # The later forecast will have narrowed, and it should still include the
    # peak size and time in its 95% credible intervals.
    cints_n2 = forecast_n2['peak_cints']
    ci_n2 = cints_n2[cints_n2['prob'] == 95]
    ci_n2_size_lower = ci_n2['sizemin'].item()
    ci_n2_size_upper = ci_n2['sizemax'].item()
    ci_n2_date_lower = datetime.datetime.strptime(
        ci_n2['timemin'].item().decode(),
        dt_format)
    ci_n2_date_upper = datetime.datetime.strptime(
        ci_n2['timemax'].item().decode(),
        dt_format)
    assert ci_n2_size_lower <= peak_size <= ci_n2_size_upper
    assert ci_n2_date_lower <= peak_date <= ci_n2_date_upper

    # The later forecast should have more accurate predictions of peak size.
    size_acc_n1 = forecast_n1['peak_size_acc']['acc']
    size_acc_n2 = forecast_n2['peak_size_acc']['acc']
    assert all(size_acc_n1 > 0.3)
    assert any(size_acc_n1 < 0.7)
    assert all(size_acc_n2 > 0.7)

    # The later forecast should have more accurate predictions of peak time.
    time_acc_n1 = forecast_n1['peak_time_acc']['acc']
    time_acc_n2 = forecast_n2['peak_time_acc']['acc']
    assert all(time_acc_n1 > 0.1)
    assert any(time_acc_n1 < 0.3)
    assert all(time_acc_n2 > 0.7)

    # Clean up: remove created files.
    os.remove(obs_file)
    os.remove(pr_file)
    os.remove(state[fs_date_n1]['forecast_file'])
    os.remove(state[fs_date_n2]['forecast_file'])


def test_seir_forecast():
    """
    Use the SEIR forecasting example to compare peak size and time predictions
    at two forecasting dates.

    Note that the observation probability is set to 0.5 (much too high) and so
    we should only obtain sensible forecasts if the observation model is able
    to use the lookup table and obtain observation probabilities from the
    ``pr-obs.ssv`` data file.
    """
    logging.basicConfig(level=logging.INFO)

    toml_file = 'seir.toml'
    obs_file = 'weekly-cases.ssv'
    pr_file = 'pr-obs.ssv'

    toml_data = pkgutil.get_data('epifx.example.seir', toml_file).decode()
    config = pypfilt.config.from_string(toml_data)

    obs_data = pkgutil.get_data('epifx.example.seir', obs_file).decode()
    with open(obs_file, mode='w') as f:
        f.write(obs_data)

    pr_data = pkgutil.get_data('epifx.example.seir', pr_file).decode()
    with open(pr_file, mode='w') as f:
        f.write(pr_data)

    forecast_from = datetime.datetime(2014, 4, 1)

    # Check that there is only one set of forecasts (i.e., only one location
    # and only one set of observation model parameters).
    forecasts = pypfilt.sweep.forecasts(config)
    forecasts = list(forecasts)
    assert len(forecasts) == 1

    # Check that forecasts were run for two forecasting dates.
    forecast = forecasts[0]
    forecast_dates = two_forecast_dates(forecast.all_observations,
                                        forecast_from)
    state = fs.run(forecast, forecast_dates)
    fs_dates = list(state.keys())
    assert len(fs_dates) == 2

    fs_date_n1 = fs_dates[0]
    fs_date_n2 = fs_dates[1]

    # Retrieve the list of observations
    obs = state[fs_date_n1]['obs']
    peak_size = max(o['value'] for o in obs)
    peak_date = [o['date'] for o in obs if o['value'] == peak_size][0]

    # Check that the peak size and date is as expected.
    assert peak_size == 2678
    assert peak_date == datetime.datetime(2014, 9, 14)

    # Compare the forecast predictions to the observed peak size and date.
    forecast_n1 = state[fs_date_n1][fs_date_n1]['summary']
    forecast_n2 = state[fs_date_n2][fs_date_n2]['summary']
    dt_format = '%Y-%m-%d %H:%M:%S'

    # The earlier forecast should include the peak size and time in its 95%
    # credible intervals.
    cints_n1 = forecast_n1['peak_cints']
    ci_n1 = cints_n1[cints_n1['prob'] == 95]
    ci_n1_size_lower = ci_n1['sizemin'].item()
    ci_n1_size_upper = ci_n1['sizemax'].item()
    ci_n1_date_lower = datetime.datetime.strptime(
        ci_n1['timemin'].item().decode(),
        dt_format)
    ci_n1_date_upper = datetime.datetime.strptime(
        ci_n1['timemax'].item().decode(),
        dt_format)
    assert ci_n1_size_lower <= peak_size <= ci_n1_size_upper
    assert ci_n1_date_lower <= peak_date <= ci_n1_date_upper

    # The later forecast will have narrowed, and it should still include the
    # peak size and time in its 95% credible intervals.
    cints_n2 = forecast_n2['peak_cints']
    ci_n2 = cints_n2[cints_n2['prob'] == 95]
    ci_n2_size_lower = ci_n2['sizemin'].item()
    ci_n2_size_upper = ci_n2['sizemax'].item()
    ci_n2_date_lower = datetime.datetime.strptime(
        ci_n2['timemin'].item().decode(),
        dt_format)
    ci_n2_date_upper = datetime.datetime.strptime(
        ci_n2['timemax'].item().decode(),
        dt_format)
    assert ci_n2_size_lower <= peak_size <= ci_n2_size_upper
    assert ci_n2_date_lower <= peak_date <= ci_n2_date_upper

    # The later forecast should have more accurate predictions of peak size.
    size_acc_n1 = forecast_n1['peak_size_acc']['acc']
    size_acc_n2 = forecast_n2['peak_size_acc']['acc']
    assert all(size_acc_n1 > 0.3)
    assert any(size_acc_n1 < 0.7)
    assert all(size_acc_n2 > 0.7)

    # The later forecast should have more accurate predictions of peak time.
    time_acc_n1 = forecast_n1['peak_time_acc']['acc']
    time_acc_n2 = forecast_n2['peak_time_acc']['acc']
    assert all(time_acc_n1 > 0.1)
    assert any(time_acc_n1 < 0.3)
    assert all(time_acc_n2 > 0.7)

    # Clean up: remove created files.
    os.remove(obs_file)
    os.remove(pr_file)
    os.remove(state[fs_date_n1]['forecast_file'])
    os.remove(state[fs_date_n2]['forecast_file'])


def test_seeiir_scalar_forecast():
    """
    Use the SEEIIR forecasting example to compare peak size and time
    predictions at two forecasting dates.

    Note that the observation probability is set to 0.5 (much too high) and so
    we should only obtain sensible forecasts if the observation model is able
    to use the lookup table and obtain observation probabilities from the
    ``pr-obs.ssv`` data file.
    """
    logging.basicConfig(level=logging.INFO)

    toml_file = 'seeiir_scalar.toml'
    obs_file = 'weekly-cases-scalar.ssv'
    pr_file = 'pr-obs-scalar.ssv'

    toml_data = pkgutil.get_data('epifx.example.seir', toml_file).decode()
    config = pypfilt.config.from_string(toml_data)

    obs_data = pkgutil.get_data('epifx.example.seir', obs_file).decode()
    with open(obs_file, mode='w') as f:
        f.write(obs_data)

    pr_data = pkgutil.get_data('epifx.example.seir', pr_file).decode()
    with open(pr_file, mode='w') as f:
        f.write(pr_data)

    forecast_from = 91

    # Check that there is only one set of forecasts (i.e., only one location
    # and only one set of observation model parameters).
    forecasts = pypfilt.sweep.forecasts(config)
    forecasts = list(forecasts)
    assert len(forecasts) == 1

    # Check that forecasts were run for two forecasting dates.
    forecast = forecasts[0]
    forecast_dates = two_forecast_times(forecast.all_observations,
                                        forecast_from)
    state = fs.run(forecast, forecast_dates)
    fs_dates = list(state.keys())
    assert len(fs_dates) == 2

    fs_date_n1 = fs_dates[0]
    fs_date_n2 = fs_dates[1]

    # Retrieve the list of observations
    obs = state[fs_date_n1]['obs']
    peak_size = max(o['value'] for o in obs)
    peak_date = [o['date'] for o in obs if o['value'] == peak_size][0]

    # Check that the peak size and date is as expected.
    assert peak_size == 2678
    assert peak_date == 257

    # Compare the forecast predictions to the observed peak size and date.
    forecast_n1 = state[fs_date_n1][fs_date_n1]['summary']
    forecast_n2 = state[fs_date_n2][fs_date_n2]['summary']

    # The earlier forecast should include the peak size and time in its 95%
    # credible intervals.
    cints_n1 = forecast_n1['peak_cints']
    ci_n1 = cints_n1[cints_n1['prob'] == 95]
    ci_n1_size_lower = ci_n1['sizemin'].item()
    ci_n1_size_upper = ci_n1['sizemax'].item()
    ci_n1_date_lower = ci_n1['timemin'].item()
    ci_n1_date_upper = ci_n1['timemax'].item()
    assert ci_n1_size_lower <= peak_size <= ci_n1_size_upper
    assert ci_n1_date_lower <= peak_date <= ci_n1_date_upper

    # The later forecast will have narrowed, and it should still include the
    # peak size and time in its 95% credible intervals.
    cints_n2 = forecast_n2['peak_cints']
    ci_n2 = cints_n2[cints_n2['prob'] == 95]
    ci_n2_size_lower = ci_n2['sizemin'].item()
    ci_n2_size_upper = ci_n2['sizemax'].item()
    ci_n2_date_lower = ci_n2['timemin'].item()
    ci_n2_date_upper = ci_n2['timemax'].item()
    assert ci_n2_size_lower <= peak_size <= ci_n2_size_upper
    assert ci_n2_date_lower <= peak_date <= ci_n2_date_upper

    # The later forecast should have more accurate predictions of peak size.
    size_acc_n1 = forecast_n1['peak_size_acc']['acc']
    size_acc_n2 = forecast_n2['peak_size_acc']['acc']
    assert all(size_acc_n1 > 0.3)
    assert any(size_acc_n1 < 0.7)
    assert all(size_acc_n2 > 0.7)

    # The later forecast should have more accurate predictions of peak time.
    time_acc_n1 = forecast_n1['peak_time_acc']['acc']
    time_acc_n2 = forecast_n2['peak_time_acc']['acc']
    assert all(time_acc_n1 > 0.1)
    assert any(time_acc_n1 < 0.3)
    assert all(time_acc_n2 > 0.7)

    # Clean up: remove created files.
    os.remove(obs_file)
    os.remove(pr_file)
    os.remove(state[fs_date_n1]['forecast_file'])
    os.remove(state[fs_date_n2]['forecast_file'])
