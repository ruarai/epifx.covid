"Test the pypfilt.cache module."

import datetime
import numpy as np
import os
import h5py
import pypfilt
import pypfilt.cache
import pypfilt.context
import pypfilt.examples
import pypfilt.examples.predation


import test_predation


def test_cache_scalar():
    cache_file = 'test_cache_scalar.hdf5'
    time_scale = pypfilt.Scalar()
    t0 = 0.0
    t1 = 15.0

    # Simulate observations from a single particle.
    sim_obs = test_predation.simulate_from_model(time_scale, t0, t1,
                                                 px_count=1)

    # Extract the observations from the summary table.
    x_tbl = sim_obs[sim_obs['unit'] == 'x'][['date', 'value']]
    y_tbl = sim_obs[sim_obs['unit'] == 'y'][['date', 'value']]
    obs = [{
        'date': time_scale.from_dtype(row['date']),
        'period': 0,
        'unit': row['unit'],
        'value': row['value'],
        'source': 'test_cache',
    } for row in sim_obs]

    # Run several forecasts to populate the cache, and then check that the
    # cache retrieves the correct state for various contexts.
    # The forecasts are, in chronological order:
    #
    # A: with original obs table and lookup table
    # B: change obs with date between A and B
    # C: original obs and lookup table
    # D: original obs and different lookup table

    params = test_predation.make_params(time_scale)
    params['time']['start'] = t0
    params['time']['until'] = t1
    params['hist']['cache_file'] = cache_file
    params['data']['obs']['x'] = x_tbl
    params['data']['obs']['y'] = y_tbl
    params['data']['lookup']['test'] = x_tbl.copy()

    def define_summary(params, obs):
        summary = pypfilt.summary.HDF5(params, obs)
        params['component']['summary'] = summary
        params['component']['summary_table'] = {
            'model_cints': pypfilt.summary.ModelCIs(probs=[0, 50, 95]),
            'obs': pypfilt.summary.Obs(),
        }

    # Forecast A: with original observation and lookup tables
    fs_time_A = t0 + 3
    define_summary(params, obs)
    pypfilt.forecast(params, [obs], [fs_time_A], filename=None)
    assert os.path.isfile(cache_file)

    # Forecast B: change obs with date between A and B
    fs_time_B = t0 + 5
    change_date = t0 + 4
    y_new = y_tbl.copy()
    orig_value = y_new['value'][y_new['date'] == change_date]
    y_new['value'][y_new['date'] == change_date] = orig_value + 0.1
    new_obs = [o.copy() for o in obs]
    for o in new_obs:
        if o['date'] == change_date and o['unit'] == 'y':
            o['value'] = o['value'] + 0.1
    params['data']['obs']['x'] = x_tbl
    params['data']['obs']['y'] = y_new
    define_summary(params, new_obs)
    pypfilt.forecast(params, [obs], [fs_time_B], filename=None)

    # Forecast C: original obs and lookup table
    fs_time_C = t0 + 7
    params['data']['obs']['x'] = x_tbl
    params['data']['obs']['y'] = y_tbl
    define_summary(params, obs)
    pypfilt.forecast(params, [obs], [fs_time_C], filename=None)

    # Forecast D: original obs and different lookup table
    fs_time_D = t0 + 9
    params['data']['lookup']['test'] = y_tbl.copy()
    define_summary(params, obs)
    pypfilt.forecast(params, [obs], [fs_time_D], filename=None)

    # Create a context object that is consistent with forecasts A and C.
    params_AC = test_predation.make_params(pypfilt.Scalar())
    params_AC['component']['summary'] = pypfilt.summary.HDF5(params, obs)
    params_AC['hist']['cache_file'] = cache_file
    params_AC['data']['obs']['x'] = x_tbl
    params_AC['data']['obs']['y'] = y_tbl
    params_AC['data']['lookup']['test'] = x_tbl.copy()
    params_AC['time'] = {
        'start': t0,
        'until': t1,
    }
    ctx_AC = pypfilt.context.Context(params_AC)

    # Create a context object that is consistent with forecast B.
    params_B = test_predation.make_params(pypfilt.Scalar())
    params_B['component']['summary'] = pypfilt.summary.HDF5(params, obs)
    params_B['hist']['cache_file'] = cache_file
    params_B['data']['obs']['x'] = x_tbl
    params_B['data']['obs']['y'] = y_new
    params_B['data']['lookup']['test'] = x_tbl.copy()
    params_B['time'] = {
        'start': t0,
        'until': t1,
    }
    ctx_B = pypfilt.context.Context(params_B)

    # Create a context object that is consistent with forecast D.
    params_D = test_predation.make_params(pypfilt.Scalar())
    params_D['component']['summary'] = pypfilt.summary.HDF5(params, obs)
    params_D['hist']['cache_file'] = cache_file
    params_D['data']['obs']['x'] = x_tbl
    params_D['data']['obs']['y'] = y_tbl
    params_D['data']['lookup']['test'] = y_tbl.copy()
    params_D['time'] = {
        'start': t0,
        'until': t1,
    }
    ctx_D = pypfilt.context.Context(params_D)

    # Check the cache results when using the same observations and lookup
    # tables as forecasts A and C.
    result = pypfilt.cache.load_state(cache_file, ctx_AC, [fs_time_A])
    assert result is not None
    assert result['start'] == fs_time_A

    result = pypfilt.cache.load_state(cache_file, ctx_AC, [fs_time_B])
    assert result is not None
    assert result['start'] == fs_time_A

    result = pypfilt.cache.load_state(cache_file, ctx_AC, [fs_time_C])
    assert result is not None
    assert result['start'] == fs_time_C

    result = pypfilt.cache.load_state(cache_file, ctx_AC, [fs_time_D])
    assert result is not None
    assert result['start'] == fs_time_C

    result = pypfilt.cache.load_state(cache_file, ctx_AC, [t0])
    assert result is None

    result = pypfilt.cache.load_state(cache_file, ctx_AC, [t1])
    assert result is not None
    assert result['start'] == fs_time_C

    # Check the cache results when using the same observations and lookup
    # tables as forecast B.
    result = pypfilt.cache.load_state(cache_file, ctx_B, [fs_time_A])
    assert result is not None
    assert result['start'] == fs_time_A

    result = pypfilt.cache.load_state(cache_file, ctx_B, [fs_time_B])
    assert result is not None
    assert result['start'] == fs_time_B

    result = pypfilt.cache.load_state(cache_file, ctx_B, [fs_time_C])
    assert result is not None
    assert result['start'] == fs_time_B

    result = pypfilt.cache.load_state(cache_file, ctx_B, [fs_time_D])
    assert result is not None
    assert result['start'] == fs_time_B

    result = pypfilt.cache.load_state(cache_file, ctx_B, [t0])
    assert result is None

    result = pypfilt.cache.load_state(cache_file, ctx_B, [t1])
    assert result is not None
    assert result['start'] == fs_time_B

    # Check the cache results when using the same observations and lookup
    # tables as forecast D.
    result = pypfilt.cache.load_state(cache_file, ctx_D, [fs_time_A])
    assert result is None

    result = pypfilt.cache.load_state(cache_file, ctx_D, [fs_time_B])
    assert result is None

    result = pypfilt.cache.load_state(cache_file, ctx_D, [fs_time_C])
    assert result is None

    result = pypfilt.cache.load_state(cache_file, ctx_D, [fs_time_D])
    assert result is not None
    assert result['start'] == fs_time_D

    result = pypfilt.cache.load_state(cache_file, ctx_D, [t0])
    assert result is None

    result = pypfilt.cache.load_state(cache_file, ctx_D, [t1])
    assert result is not None
    assert result['start'] == fs_time_D

    os.remove(cache_file)


def test_cache_datetime():
    cache_file = 'test_cache_datetime.hdf5'
    time_scale = pypfilt.Datetime()
    t0 = datetime.datetime(2017, 5, 1)
    t1 = t0 + datetime.timedelta(days=15)

    # Simulate observations from a single particle.
    sim_obs = test_predation.simulate_from_model(time_scale, t0, t1,
                                                 px_count=1)

    # Convert the date column from strings to datetime values.
    rows = [(time_scale.from_dtype(row['date']),
             row['unit'], row['value'])
            for row in sim_obs]
    dtype = [('date', time_scale.native_dtype()),
             ('unit', h5py.string_dtype(encoding='utf-8')),
             ('value', np.float_)]
    obs_tbl = np.array(rows, dtype=dtype)

    # Extract the observations from the summary table.
    x_tbl = obs_tbl[obs_tbl['unit'] == 'x'][['date', 'unit', 'value']]
    y_tbl = obs_tbl[obs_tbl['unit'] == 'y'][['date', 'unit', 'value']]
    obs = [{
        'date': row['date'],
        'period': 0,
        'unit': row['unit'],
        'value': row['value'],
        'source': 'test_cache',
    } for row in obs_tbl]

    sim_obs = None
    obs_tbl = None

    params = test_predation.make_params(time_scale)
    params['time']['start'] = t0
    params['time']['until'] = t1
    params['hist']['cache_file'] = cache_file
    params['data']['obs']['x'] = x_tbl
    params['data']['obs']['y'] = y_tbl
    params['data']['lookup']['test'] = x_tbl.copy()

    def define_summary(params, obs):
        summary = pypfilt.summary.HDF5(params, obs)
        params['component']['summary'] = summary
        params['component']['summary_table'] = {
            'model_cints': pypfilt.summary.ModelCIs(probs=[0, 50, 95]),
            'obs': pypfilt.summary.Obs(),
        }

    # Forecast A: with original observation and lookup tables
    fs_time_A = t0 + datetime.timedelta(days=3)
    define_summary(params, obs)
    pypfilt.forecast(params, [obs], [fs_time_A], filename=None)
    assert os.path.isfile(cache_file)

    # Forecast B: change obs with date between A and B
    fs_time_B = t0 + datetime.timedelta(days=5)
    change_date = t0 + datetime.timedelta(days=4)
    y_new = y_tbl.copy()
    orig_value = y_new['value'][y_new['date'] == change_date]
    y_new['value'][y_new['date'] == change_date] = orig_value + 0.1
    new_obs = [o.copy() for o in obs]
    for o in new_obs:
        if o['date'] == change_date and o['unit'] == 'y':
            o['value'] = o['value'] + 0.1
    params['data']['obs']['x'] = x_tbl
    params['data']['obs']['y'] = y_new
    define_summary(params, new_obs)
    pypfilt.forecast(params, [obs], [fs_time_B], filename=None)

    # Forecast C: original obs and lookup table
    fs_time_C = t0 + datetime.timedelta(days=7)
    params['data']['obs']['x'] = x_tbl
    params['data']['obs']['y'] = y_tbl
    define_summary(params, obs)
    pypfilt.forecast(params, [obs], [fs_time_C], filename=None)

    # Forecast D: original obs and different lookup table
    fs_time_D = t0 + datetime.timedelta(days=9)
    params['data']['lookup']['test'] = y_tbl.copy()
    define_summary(params, obs)
    pypfilt.forecast(params, [obs], [fs_time_D], filename=None)

    # Create a context object that is consistent with forecasts A and C.
    params_AC = test_predation.make_params(pypfilt.Datetime())
    params_AC['component']['summary'] = pypfilt.summary.HDF5(params, obs)
    params_AC['hist']['cache_file'] = cache_file
    params_AC['data']['obs']['x'] = x_tbl
    params_AC['data']['obs']['y'] = y_tbl
    params_AC['data']['lookup']['test'] = x_tbl.copy()
    params_AC['time'] = {
        'start': t0,
        'until': t1,
    }
    ctx_AC = pypfilt.context.Context(params_AC)

    # Create a context object that is consistent with forecast B.
    params_B = test_predation.make_params(pypfilt.Datetime())
    params_B['component']['summary'] = pypfilt.summary.HDF5(params, obs)
    params_B['hist']['cache_file'] = cache_file
    params_B['data']['obs']['x'] = x_tbl
    params_B['data']['obs']['y'] = y_new
    params_B['data']['lookup']['test'] = x_tbl.copy()
    params_B['time'] = {
        'start': t0,
        'until': t1,
    }
    ctx_B = pypfilt.context.Context(params_B)

    # Create a context object that is consistent with forecast D.
    params_D = test_predation.make_params(pypfilt.Datetime())
    params_D['component']['summary'] = pypfilt.summary.HDF5(params, obs)
    params_D['hist']['cache_file'] = cache_file
    params_D['data']['obs']['x'] = x_tbl
    params_D['data']['obs']['y'] = y_tbl
    params_D['data']['lookup']['test'] = y_tbl.copy()
    params_D['time'] = {
        'start': t0,
        'until': t1,
    }
    ctx_D = pypfilt.context.Context(params_D)

    # Check the cache results when using the same observations and lookup
    # tables as forecasts A and C.
    result = pypfilt.cache.load_state(cache_file, ctx_AC, [fs_time_A])
    assert result is not None
    assert result['start'] == fs_time_A

    result = pypfilt.cache.load_state(cache_file, ctx_AC, [fs_time_B])
    assert result is not None
    assert result['start'] == fs_time_A

    result = pypfilt.cache.load_state(cache_file, ctx_AC, [fs_time_C])
    assert result is not None
    assert result['start'] == fs_time_C

    result = pypfilt.cache.load_state(cache_file, ctx_AC, [fs_time_D])
    assert result is not None
    assert result['start'] == fs_time_C

    result = pypfilt.cache.load_state(cache_file, ctx_AC, [t0])
    assert result is None

    result = pypfilt.cache.load_state(cache_file, ctx_AC, [t1])
    assert result is not None
    assert result['start'] == fs_time_C

    # Check the cache results when using the same observations and lookup
    # tables as forecast B.
    result = pypfilt.cache.load_state(cache_file, ctx_B, [fs_time_A])
    assert result is not None
    assert result['start'] == fs_time_A

    result = pypfilt.cache.load_state(cache_file, ctx_B, [fs_time_B])
    assert result is not None
    assert result['start'] == fs_time_B

    result = pypfilt.cache.load_state(cache_file, ctx_B, [fs_time_C])
    assert result is not None
    assert result['start'] == fs_time_B

    result = pypfilt.cache.load_state(cache_file, ctx_B, [fs_time_D])
    assert result is not None
    assert result['start'] == fs_time_B

    result = pypfilt.cache.load_state(cache_file, ctx_B, [t0])
    assert result is None

    result = pypfilt.cache.load_state(cache_file, ctx_B, [t1])
    assert result is not None
    assert result['start'] == fs_time_B

    # Check the cache results when using the same observations and lookup
    # tables as forecast D.
    result = pypfilt.cache.load_state(cache_file, ctx_D, [fs_time_A])
    assert result is None

    result = pypfilt.cache.load_state(cache_file, ctx_D, [fs_time_B])
    assert result is None

    result = pypfilt.cache.load_state(cache_file, ctx_D, [fs_time_C])
    assert result is None

    result = pypfilt.cache.load_state(cache_file, ctx_D, [fs_time_D])
    assert result is not None
    assert result['start'] == fs_time_D

    result = pypfilt.cache.load_state(cache_file, ctx_D, [t0])
    assert result is None

    result = pypfilt.cache.load_state(cache_file, ctx_D, [t1])
    assert result is not None
    assert result['start'] == fs_time_D

    os.remove(cache_file)
