"""Test the pypfilt.examples.predation module."""

import datetime
import numpy as np
import os
import os.path
import pypfilt
import pypfilt.examples
import pypfilt.examples.predation
import pypfilt.sweep
import pytest


def make_params(time_scale, max_days=14):
    """Define the default simulation parameters for this model."""
    px_count = 1000
    seed = 42
    obs_sdev = 0.2
    params = pypfilt.examples.predation.make_params(px_count, seed, obs_sdev)
    params['component']['time'] = time_scale

    # Record credible intervals for model parameters and state variables.
    probs = [0, 50, 95]
    params['summary'].setdefault('tables', {}).setdefault('model_cints', {})
    params['summary']['tables']['model_cints']['credible_intervals'] = probs

    return params


def fs_scalar_time():
    # Define the simulation period and forecasting times.
    params = make_params(pypfilt.Scalar())
    t0 = 0.0
    t1 = 15.0
    fs_times = [1.0, 3.0, 5.0, 7.0, 9.0]

    # Generate noisy observations.
    params['time']['start'] = t0
    params['time']['until'] = t1
    obs = pypfilt.examples.predation.make_observations(params)

    # Define the summary tables to be saved to disk.
    summary = pypfilt.summary.HDF5(params, obs)
    params['component']['summary'] = summary
    params['component']['summary_table'] = {
        'model_cints': pypfilt.summary.ModelCIs(),
        'obs': pypfilt.summary.Obs(),
    }

    # Run the forecast simulations.
    params['time']['start'] = t0
    params['time']['until'] = t1
    return pypfilt.forecast(params, [obs], fs_times, filename=None)


def test_fs_scale_time_summarise_estimation():
    """
    This test triggers the summarise bug that was fixed by commit 1c4af76.

    This required modifying how pypfilt.cache.__data_match compares data keys,
    because no lookup tables are defined.
    """
    # Define the simulation period and forecasting times.
    cache_file = 'test_fs_scale_time_summarise_estimation.hdf5'
    if os.path.isfile(cache_file):
        os.remove(cache_file)

    params = make_params(pypfilt.Scalar(), max_days=3)
    params['summary']['only_forecasts'] = False
    params['hist']['cache_file'] = cache_file
    t0 = 0.0
    t1 = 15.0
    fs_times = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0]

    # Generate noisy observations.
    params['time']['start'] = t0
    params['time']['until'] = t1
    obs = pypfilt.examples.predation.make_observations(params)

    # Define the summary tables to be saved to disk.
    summary = pypfilt.summary.HDF5(params, obs)
    params['component']['summary'] = summary
    params['component']['summary_table'] = {
        'model_cints': pypfilt.summary.ModelCIs(),
        'obs': pypfilt.summary.Obs(),
    }

    # Run the forecast simulations.
    params['time']['start'] = t0
    params['time']['until'] = t1
    for fs_time in fs_times:
        pypfilt.forecast(params, [obs], [fs_time], filename=None)

    assert os.path.isfile(cache_file)
    os.remove(cache_file)


def fs_date_time():
    # Define the simulation period and forecasting times.
    params = make_params(pypfilt.Datetime())
    t0 = datetime.datetime(2017, 5, 1)
    t1 = t0 + datetime.timedelta(days=15)
    fs_times = [t0 + datetime.timedelta(days=t)
                for t in [1.0, 3.0, 5.0, 7.0, 9.0]]

    # Generate noisy observations.
    params['time']['start'] = t0
    params['time']['until'] = t1
    obs = pypfilt.examples.predation.make_observations(params)

    # Define the summary tables to be saved to disk.
    summary = pypfilt.summary.HDF5(params, obs)
    params['component']['summary'] = summary
    params['component']['summary_table'] = {
        'model_cints': pypfilt.summary.ModelCIs(),
        'obs': pypfilt.summary.Obs(),
    }

    # Run the forecast simulations.
    params['time']['start'] = t0
    params['time']['until'] = t1
    return pypfilt.forecast(params, [obs], fs_times, filename=None)


def cmp_summaries(st_sum, dt_sum):
    assert set(st_sum) == set(dt_sum)
    for name, st_tbl in st_sum.items():
        dt_tbl = dt_sum[name]
        st_cols = st_tbl.dtype.descr
        dt_cols = dt_tbl.dtype.descr
        for ix, st_col in enumerate(st_cols):
            dt_col = dt_cols[ix]
            assert st_col[0] == dt_col[0]
            if st_col[1] == dt_col[1]:
                assert (st_tbl[st_col[0]] == dt_tbl[dt_col[0]]).all()
    return True


def test_time_scales():
    st_tbls = fs_scalar_time()
    dt_tbls = fs_date_time()
    # Compare observations.
    for col in ['unit', 'value', 'source', 'period']:
        st_vals = [o[col] for o in st_tbls['obs']]
        dt_vals = [o[col] for o in dt_tbls['obs']]
        assert st_vals == dt_vals
    # Compare estimation run summaries.
    st_sum = st_tbls['complete']['summary']
    dt_sum = dt_tbls['complete']['summary']
    assert cmp_summaries(st_sum, dt_sum)
    # Compare summaries for each forecasting date.
    dt_fs_dates = [k for k in dt_tbls if isinstance(k, datetime.datetime)]
    dt_t0 = datetime.datetime(2017, 5, 1)
    for dt_fs in sorted(dt_fs_dates):
        st_fs = (dt_fs - dt_t0).days
        assert st_fs in st_tbls
        st_sum = st_tbls[st_fs]['summary']
        dt_sum = dt_tbls[dt_fs]['summary']
        assert cmp_summaries(st_sum, dt_sum)


def test_main():
    pypfilt.examples.predation.main()
    output_files = ['predation.hdf5', 'predation_forecasts.png',
                    'predation_params.png']
    for filename in output_files:
        assert os.path.isfile(filename)
        os.remove(filename)


@pytest.fixture(params=[True, False])
def use_cache_file(request):
    return request.param


def test_scalar_time_no_window(use_cache_file):
    """
    Run forecasts where the state history matrix contains the entire
    simulation period, to ensure that the summary outputs cover the entire
    simulation period too.
    """
    # Define the simulation period and forecasting times.
    params = make_params(pypfilt.Scalar())
    t0 = 0.0
    t1 = 7.0
    fs_times = [3.0]

    if use_cache_file:
        cache_file = 'test_scalar_time_no_window.hdf5'
        params['hist']['cache_file'] = cache_file
        if os.path.isfile(cache_file):
            os.remove(cache_file)

    # Generate noisy observations.
    params['time']['start'] = t0
    params['time']['until'] = t1
    obs = pypfilt.examples.predation.make_observations(params)

    # Define the summary tables to be saved to disk.
    summary = pypfilt.summary.HDF5(params, obs)
    params['component']['summary'] = summary
    params['component']['summary_table'] = {
        'model_cints': pypfilt.summary.ModelCIs(),
        'obs': pypfilt.summary.Obs(),
    }

    # Run the forecast simulations twice; the first time will create the cache
    # and the second time will use the cache (if we're using a cache file).
    params['time']['start'] = t0
    params['time']['until'] = t1
    state0 = pypfilt.forecast(params, [obs], fs_times, filename=None)
    state1 = pypfilt.forecast(params, [obs], fs_times, filename=None)

    # Check that the first forecasts required an estimation run, and that this
    # is where the observations were recorded.
    assert 'complete' in state0
    assert 'obs' in state0['complete']['summary']
    assert 'obs' not in state0[3.0]['summary']

    # Check that the model CIs were recorded in each simulation.
    assert 'model_cints' in state0['complete']['summary']
    assert 'model_cints' in state0[3.0]['summary']
    assert 'model_cints' in state1[3.0]['summary']

    # Check that the model CIs were recorded as expected.
    assert (
        set(np.unique(state0['complete']['summary']['model_cints']['prob']))
        == {0, 50, 95})

    if use_cache_file:
        # Check that the second forecasts did not require an estimation run,
        # and that the observations were recorded.
        assert 'complete' not in state1
        assert 'obs' in state1[3.0]['summary']
    else:
        # Check that the first forecasts required an estimation run, and that
        # this is where the observations were recorded.
        assert 'complete' in state1
        assert 'obs' in state1['complete']['summary']
        assert 'obs' not in state1[3.0]['summary']
        assert 'model_cints' in state1['complete']['summary']

    def check_table_date_ranges(state, fs_key, t_start):
        for table_name in state[fs_key]['summary']:
            table = state[fs_key]['summary'][table_name]

            # Determine the appropriate date range for each table.
            if table_name == 'obs':
                # The first observation occurs at time t = 1.0, with the
                # exception of simulated observations, which start at t0.
                min_date = t0
            else:
                # The model credible intervals begin at the starting time
                # (this is 0.0 for the estimation run, and the forecasting
                # time for forecasts).
                min_date = t_start

            # The estimation run will stop at the final forecasting date, but
            # all of the observations (i.e., up to the end of the simulation
            # period) will still be recorded.
            if fs_key == 'complete' and table_name != 'obs':
                max_date = max(fs_times)
            else:
                max_date = t1

            assert min(table['date']) >= min_date
            assert max(table['date']) == max_date

    check_table_date_ranges(state0, 'complete', t0)
    check_table_date_ranges(state0, 3.0, 3.0)
    if not use_cache_file:
        check_table_date_ranges(state1, 'complete', t0)
    check_table_date_ranges(state1, 3.0, 3.0)

    if use_cache_file:
        assert os.path.isfile(cache_file)
        os.remove(cache_file)


def test_config():
    obs_x_file = 'predation-counts-x.ssv'
    obs_y_file = 'predation-counts-y.ssv'

    toml_data = pypfilt.examples.predation.example_toml_data()
    config = pypfilt.config.from_string(toml_data)

    obs_x_data = pypfilt.examples.predation.example_obs_x_data()
    with open(obs_x_file, 'w') as f:
        f.write(obs_x_data)

    obs_y_data = pypfilt.examples.predation.example_obs_y_data()
    with open(obs_y_file, 'w') as f:
        f.write(obs_y_data)

    forecasts = pypfilt.sweep.forecasts(config)
    forecasts = list(forecasts)
    assert len(forecasts) == 1

    forecast = forecasts[0]
    forecast_times = [1.0, 3.0, 5.0, 7.0, 9.0]

    state = pypfilt.forecast(forecast.params,
                             forecast.observation_streams,
                             forecast_times,
                             filename=None)

    assert 'complete' in state
    assert 'obs' in state
    for t in forecast_times:
        assert t in state
        # Compare the credible intervals for x(t) and y(t) to the forecast
        # credible intervals; the forecast CIs should be wider because they
        # include the observation variance.
        assert 'summary' in state[t]
        tables = state[t]['summary']
        assert 'model_cints' in tables
        assert 'forecasts' in tables
        model_cints = tables['model_cints']
        forecasts = tables['forecasts']
        probs = np.unique(forecasts['prob'])
        probs = [pr for pr in probs if pr > 0]
        for name in np.unique(forecasts['unit']):
            for pr in probs:
                m_ci = model_cints[
                    np.logical_and(model_cints['name'] == name,
                                   model_cints['prob'] == pr)]
                m_fs = forecasts[
                    np.logical_and(forecasts['unit'] == name,
                                   forecasts['prob'] == pr)]
                assert m_ci.shape == m_fs.shape
                assert np.all(m_ci['ymin'] > m_fs['ymin'])
                assert np.all(m_ci['ymax'] < m_fs['ymax'])

    for data_file in [obs_x_file, obs_y_file]:
        assert os.path.isfile(data_file)
        os.remove(data_file)


def test_simulated_obs():
    # Define the simulation period and forecasting times.
    params = make_params(pypfilt.Scalar())
    t0 = 0.0
    t1 = 15.0
    fs_times = [1.0, 3.0, 5.0, 7.0, 9.0]

    # Generate noisy observations.
    params['time']['start'] = t0
    params['time']['until'] = t1
    obs = pypfilt.examples.predation.make_observations(params)

    # Define the summary tables to be saved to disk.
    summary = pypfilt.summary.HDF5(params, obs)
    params['component']['summary'] = summary
    params['component']['summary_table'] = {
        'model_cints': pypfilt.summary.ModelCIs(),
        'obs': pypfilt.summary.Obs(),
        'simulated_obs': pypfilt.summary.SimulatedObs()
    }

    # Run the forecast simulations.
    params['time']['start'] = t0
    params['time']['until'] = t1
    state = pypfilt.forecast(params, [obs], fs_times, filename=None)

    # Check the simulated observations for each forecast.
    for fs_time in fs_times:
        assert fs_time in state
        tables = state[fs_time]['summary']
        assert 'simulated_obs' in tables
        sim_obs = tables['simulated_obs']
        assert all(sim_obs['fs_date'] == fs_time)
        assert all(sim_obs['date'] >= fs_time)
        num_particles = np.sum(sim_obs['date'] == fs_time)
        assert len(sim_obs) == num_particles * int(t1 - fs_time + 1)


def simulate_from_model(time_scale, t0, t1, px_count):
    seed = 42
    obs_sdev = 0.2

    params = pypfilt.examples.predation.make_params(px_count, seed, obs_sdev)
    params['component']['time'] = time_scale
    params['time']['start'] = t0
    params['time']['until'] = t1

    # Define the ground truth and construct the corresponding priors.
    x0 = 0.9
    y0 = 0.25
    alpha = 2/3
    beta = 4/3
    gamma = 1
    delta = 1
    params['model']['prior'] = {
        'x': lambda r, size=None: x0 * np.ones(size),
        'y': lambda r, size=None: y0 * np.ones(size),
        'alpha': lambda r, size=None: alpha * np.ones(size),
        'beta': lambda r, size=None: beta * np.ones(size),
        'gamma': lambda r, size=None: gamma * np.ones(size),
        'delta': lambda r, size=None: delta * np.ones(size),
    }

    # Simulate the observations from this model.
    return pypfilt.simulate_from_model(params)


def test_simulate_from_model():
    # Simulate the observations from this model.
    t0 = 0.0
    t1 = 15.0
    px_count = 1
    sim_obs = simulate_from_model(pypfilt.Scalar(), t0, t1, px_count)

    # Check that the output appears sensible.
    num_obs_dates = int(t1 - t0 + 1)
    assert isinstance(sim_obs, np.ndarray)
    assert len(sim_obs) == 2 * num_obs_dates
    assert all(sim_obs['fs_date'] == 0.0)
    assert sum(sim_obs['unit'] == 'x') == num_obs_dates
    assert sum(sim_obs['unit'] == 'y') == num_obs_dates
    assert all(sim_obs['value'] > -0.5)
    assert all(sim_obs['value'] < 2.5)
    for t in range(num_obs_dates):
        obs_at_t = sim_obs[sim_obs['date'] == float(t)]
        assert len(obs_at_t) == 2
        assert any(obs_at_t['unit'] == 'x')
        assert any(obs_at_t['unit'] == 'y')

    # pypfilt.examples.predation.save_scalar_observations(sim_obs)
