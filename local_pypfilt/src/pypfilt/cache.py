"""Particle filter state cache."""

import atexit
import logging
import h5py
import numpy as np
import os.path
import signal
import tempfile


def default(ctx, forecast_dates):
    """
    Return the default (i.e., empty) cached state.

    If no cache file is defined in ``ctx.params['hist']['cache_file']``, this
    will create a temporary cache file that will be automatically deleted when
    the process terminates (either normally, or in response to the SIGTERM
    signal).

    :param ctx: The simulation context.
    :param forecast_dates: The dates at which forecasts will be run.
    :returns: A dictionary with keys:

        - ``'state'``: Either a simulation state (see
          :func:`pypfilt.pfilter.run`) or ``None`` if there is no cached state
          from which to resume.
        - ``'start'``: Either the time from which to begin the simulation, or
          ``None`` if there is no cached state.
        - ``'end'``: Either the time at which to end the simulation, or
          ``None`` if there is no cached state.
        - ``'fs_dates'``: The dates at which forecasts should be generated.
        - ``'save_to'``: The filename for saving the particle history matrix.
        - ``'clean'``: A cleanup function to remove any temporary files, and
          which will have been registered to execute at termination.
    """
    cache_file = None
    if 'cache_file' in ctx.params['hist']:
        if ctx.params['hist']['cache_file'] is not None:
            cache_file = os.path.join(ctx.params['out_dir'],
                                      ctx.params['hist']['cache_file'])

    if cache_file is None:
        cache_file, clean_fn = temporary_cache(ctx)
    else:
        def clean_fn():
            pass

    # The default, should there be no suitable cached state.
    result = {'state': None, 'start': None, 'end': max(forecast_dates),
              'fs_dates': forecast_dates,
              'save_to': cache_file, 'clean': clean_fn}
    return result


def temporary_cache(ctx):
    logger = logging.getLogger(__name__)
    tmp_dir = tempfile.mkdtemp(dir=ctx.params['tmp_dir'])
    tmp_file = os.path.join(tmp_dir, "history.hdf5")

    # Ensure these files are always deleted upon *normal* termination.
    atexit.register(__cleanup, files=[tmp_file], dirs=[tmp_dir])

    # Ensure these files are always deleted when killed by SIGTERM.
    def clean_at_terminate(signal_num, stack_frame):
        __cleanup(files=[tmp_file], dirs=[tmp_dir])
        os._exit(0)

    signal.signal(signal.SIGTERM, clean_at_terminate)

    logger.debug("Temporary file for history matrix: '{}'".format(
        tmp_file))

    def clean():
        __cleanup(files=[tmp_file], dirs=[tmp_dir])

    return (tmp_file, clean)


def __cleanup(files, dirs):
    """Delete temporary files and directories.
    This is intended for use with ``atexit.register()``.

    :param files: The list of files to delete.
    :param dirs: The list of directories to delete (*after* all of the files
        have been deleted). Note that these directories must be empty in order
        to be deleted.
    """
    logger = logging.getLogger(__name__)

    for tmp_file in files:
        if os.path.isfile(tmp_file):
            try:
                os.remove(tmp_file)
                logger.debug("Deleted '{}'".format(tmp_file))
            except OSError as e:
                msg = "Can not delete '{}': {}".format(tmp_file, e.strerror)
                logger.warning(msg)
        elif os.path.exists(tmp_file):
            logger.warning("'{}' is not a file".format(tmp_file))
        else:
            logger.debug("File '{}' already deleted".format(tmp_file))

    for tmp_dir in dirs:
        if os.path.isdir(tmp_dir):
            try:
                os.rmdir(tmp_dir)
                logger.debug("Deleted '{}'".format(tmp_dir))
            except OSError as e:
                msg = "Can not delete '{}': {}".format(tmp_dir, e.strerror)
                logger.warning(msg)
        elif os.path.exists(tmp_dir):
            logger.warning("'{}' is not a directory".format(tmp_dir))
        else:
            logger.debug("Directory '{}' already deleted".format(tmp_dir))


def save_state(cache_file, ctx, when, **kwargs):
    """
    Save the particle history matrix to a cache file, to allow future
    forecasting runs to resume from this point.

    :param cache_file: The name of the cache file.
    :param ctx: The simulation context.
    :param when: The current simulation time.
    :param \\**kwargs: The data sets to store in the cached state; at a
        minimum this should include ``'hist'`` (the particle history matrix)
        and ``'offset'`` (the index of the current time-step in the particle
        history matrix).

    :Examples:

    .. code:: python

       cache_file = 'cache.hdf5'
       cache.save_state(cache_file, context, current_time,
                        offset=np.int32(hist_ix),
                        hist=np.float64(hist))
    """
    with h5py.File(cache_file, 'a') as f:
        when_str = ctx.component['time'].to_unicode(when)
        grp = f.require_group(when_str)
        # Replace any previously-saved state for this time.
        for (name, value) in kwargs.items():
            if name in grp:
                del grp[name]
            grp.create_dataset(name, data=value)
        # Save all of the data tables.
        data_grp = grp.require_group('data')
        __save_data_tables(ctx, f, data_grp, ctx.data)
        # Save the state of the summary object.
        sum_grp = grp.require_group('summary')
        ctx.component['summary'].save_state(sum_grp)


def __save_data_tables(ctx, hdf5_file, group, data, path=None):
    if path is None:
        path = []
    for (name, value) in data.items():
        if isinstance(value, dict):
            if value:
                # Nested, non-empty dictionary.
                subgroup = group.require_group(name)
                path.append(name)
                __save_data_tables(ctx, hdf5_file, subgroup, value, path)
        elif isinstance(value, np.ndarray):
            value = __convert_time_df(value, ctx.component['time'],
                                      col='date')
            if name in group:
                del group[name]
            group.create_dataset(name, data=value)
        else:
            raise ValueError('Invalid data table {}.{} has type {}'.format(
                '.'.join(path), name, type(value)))


def __is_time_df(df, col):
    return (df.dtype.names
            and col in df.dtype.names
            and df.dtype[col].type == np.object_)


def __convert_time_df(df, time, col):
    if not __is_time_df(df, col):
        return df

    descr = list(df.dtype.names)
    date_ix = None
    for ix in range(len(descr)):
        name = descr[ix]
        if name == 'date':
            descr[ix] = time.dtype('date')
            date_ix = ix
        else:
            descr[ix] = (name, df.dtype[name])

    rows = [list(row.tolist()) for row in df]
    for row in rows:
        row[date_ix] = time.to_dtype(row[date_ix])
    rows = [tuple(row) for row in rows]
    return np.array(rows, dtype=descr)


def __is_obj_df(df, time, col):
    # NOTE: use np.dtype() to convert the time dtype into an appropriate form
    # for np.issubdtype(). In a future version of NumPy, np.issubdtype() will
    # stop downcasting dtype-like arguments. For example, 'float64' was
    # translated into np.floating rather than np.float64, and so the following
    # would return True:
    # >>> issubdtype(np.float32, 'float64')
    # So we should pass actual dtype objects, rather than dtype-like values.
    return (df.dtype.names
            and col in df.dtype.names
            and np.issubdtype(df.dtype[col].type,
                              np.dtype(time.dtype('ignore')[1]))
            and np.issubdtype(np.dtype(time.dtype('ignore')[1]),
                              df.dtype[col].type))


def __convert_obj_df(df, time, col):
    if not __is_obj_df(df, time, col):
        return df

    descr = list(df.dtype.names)
    date_ix = None
    for ix in range(len(descr)):
        name = descr[ix]
        if name == 'date':
            descr[ix] = ('date', np.dtype(time.native_dtype()))
            date_ix = ix
        else:
            descr[ix] = (name, df.dtype[name])

    rows = [list(row.tolist()) for row in df]
    for row in rows:
        row[date_ix] = time.from_dtype(row[date_ix])
    rows = [tuple(row) for row in rows]
    return np.array(rows, dtype=descr)


def load_state(cache_file, ctx, forecast_dates):
    """
    Load the particle history matrix from a cache file, allowing forecasting
    runs to resume at the point of the first updated/new observation.

    :param cache_file: The name of the cache file.
    :param ctx: The simulation context.
    :param forecast_dates: The dates at which forecasts will be run.

    :returns: Either ``None``, if there was no suitable cached state, or a
        dictionary with following keys:

        - ``'start'``: The time from which to begin the simulation.
        - ``'end'``: The time from which to end the simulation.
        - ``'state'``: A dictionary that contains the following keys:

            - ``'hist'``: The particle history matrix.
            - ``'offset'``: The index of the current time-step in the particle
              history matrix.

    Note that if the cache file already contains a suitable state for each of
    the provided forecast dates, this will return a dictionary as described
    above, where the ``'start'`` and ``'end'`` values are the same (i.e.,
    there is no need to run an estimation pass).
    """
    logger = logging.getLogger(__name__)
    logger.debug('Searching for cached state in {}'.format(cache_file))

    if not os.path.exists(cache_file):
        logger.debug("Missing cache file: '{}'".format(cache_file))
        return None

    try:
        with h5py.File(cache_file, 'r') as f:
            logger.debug("Reading cache file: '{}'".format(cache_file))
            return __find_most_recent_date(ctx, forecast_dates, f)
    except IOError:
        logger.debug("Could not read cache file: '{}'".format(cache_file))
        return None


def load_state_at_time(cache_file, ctx, when):
    """
    Load the particle history matrix from a cache file at a specific
    simulation time.

    :param cache_file: The name of the cache file.
    :param ctx: The simulation context.
    :param when: The simulation time at which the particle history matrix was
        cached.

    :returns: Either ``None``, if there was no cached result for the given
        simulation time, or a dictionary as per :func:`load_state`.
    """
    logger = logging.getLogger(__name__)

    if not os.path.exists(cache_file):
        logger.debug("Missing cache file: '{}'".format(cache_file))
        return None

    time = ctx.component['time']
    summary = ctx.component['summary']
    model = ctx.component['model']

    try:
        with h5py.File(cache_file, 'r') as f:
            logger.debug("Reading cache file: '{}'".format(cache_file))
            when_str = time.to_unicode(when)
            if when_str not in f:
                logger.debug('Cache file has no entry for "{}"'
                             .format(when_str))
            group = f[when_str]
            result = {
                'start': when,
                'state': {
                    'hist': group['hist'][()],
                    'offset': group['offset'][()]
                },
            }
            summary.load_state(group['summary'])
            # Inform the model that we're resuming from a cached state.
            model.resume_from_cache(ctx)
            return result
    except IOError:
        logger.debug("Could not read cache file: '{}'".format(cache_file))
        return None


def __find_most_recent_date(ctx, forecast_dates, hdf5_file):
    time = ctx.component['time']
    summary = ctx.component['summary']
    model = ctx.component['model']

    # Starting from the earliest forecasting date, identify the forecasting
    # dates for which there is a matching cached state. This search stops when
    # it encounters the first forecasting date for which there is no matching
    # cached state, because the estimation pass will need to start from the
    # previous forecasting date.
    forecast_dates = sorted(forecast_dates)
    forecast_matches = []
    for fs_date in forecast_dates:
        group_name = time.to_unicode(fs_date)
        if group_name not in hdf5_file:
            break
        group = hdf5_file[group_name]
        if __data_match(time, ctx.data, group['data'], fs_date):
            forecast_matches.append(fs_date)
        else:
            break

    # If we have a matching cached state for the earliest forecasting date,
    # and possibly for subsequent forecasting dates, we can start the
    # estimation pass at the latest of these forecasting dates.
    #
    # If there is a matching cached state for all of the forecasting dates,
    # the 'start' and 'end' dates will be identical and the estimation pass
    # will be avoided entirely.
    if forecast_matches:
        cache_date = forecast_matches[-1]
        group_name = time.to_unicode(cache_date)
        group = hdf5_file[group_name]
        result = {
            'start': cache_date,
            'end': max(forecast_dates),
            'state': {
                'hist': group['hist'][()],
                'offset': group['offset'][()]
            },
        }
        summary.load_state(group['summary'])
        # Inform the model that we're resuming from a cached state.
        model.resume_from_cache(ctx)
        return result

    # Build a table that maps the cache date to the group name.
    cache_table = {time.from_unicode(group): group for group in hdf5_file
                   if group != 'obs'}
    # Sort cached dates from newest to oldest.
    cache_dates = reversed(sorted([date for date in cache_table.keys()
                                   if date <= min(forecast_dates)]))

    # Starting with the newest date, find the first cached state that is
    # consistent with all of the data tables in ctx.data.
    for cache_date in cache_dates:
        group_name = cache_table[cache_date]
        group = hdf5_file[group_name]
        if __data_match(time, ctx.data, group['data'], cache_date):
            result = {
                'start': cache_date,
                'end': max(forecast_dates),
                'state': {
                    'hist': group['hist'][()],
                    'offset': group['offset'][()]
                },
            }
            summary.load_state(group['summary'])
            # Inform the model that we're resuming from a cached state.
            model.resume_from_cache(ctx)
            return result

    return None


def __data_match(time, ctx_data, cache_data, cache_date, obs=False):
    # NOTE: because __save_data_tables() only creates sub-groups when it
    # encounters a *non-empty* dictionary, we need to ignore empty
    # dictionaries in the context data when comparing keys here.
    # A straight `ctx_data[k] != {}` comparison raises warnings when the value
    # is a numpy array, so we first need to check whether the value is a
    # dictionary.
    ctx_keys = set(k for k in ctx_data.keys()
                   if (not isinstance(ctx_data[k], dict))
                   or ctx_data[k] != {})
    data_keys = set(cache_data.keys())
    if ctx_keys != data_keys:
        return False

    for key in ctx_keys:
        ctx_val = ctx_data[key]
        cache_val = cache_data[key]
        if isinstance(ctx_val, dict) and isinstance(cache_val, h5py.Group):
            # Compare nested group of data tables.
            obs = (not obs) and key == 'obs'
            if not __data_match(time, ctx_val, cache_val, cache_date, obs):
                return False
        elif (isinstance(ctx_val, np.ndarray)
              and isinstance(cache_val, h5py.Dataset)):
            # Compare data tables.
            cache_tbl = cache_val[()]

            if obs:
                # NOTE: we only want to compare observations up to the cache
                # date.
                if __is_obj_df(cache_tbl, time, col='date'):
                    if cache_tbl.dtype != ctx_val.dtype:
                        cache_tbl = __convert_obj_df(cache_tbl, time,
                                                     col='date')
                    ctx_sub = ctx_val[ctx_val['date'] <= cache_date]
                    cache_sub = cache_tbl[cache_tbl['date'] <= cache_date]
                    is_match = np.array_equal(ctx_sub, cache_sub)
                else:
                    is_match = False
            else:
                ctx_tbl = __convert_time_df(ctx_val, time, col='date')
                is_match = np.array_equal(ctx_tbl, cache_tbl)

            if not is_match:
                return False
        else:
            # NOTE: mismatch between context and cache group structure.
            return False

    return True
