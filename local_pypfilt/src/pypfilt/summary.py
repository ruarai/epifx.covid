"""Calculation of simulation summary statistics."""

import abc
import h5py
import importlib
import locale
import logging
import numpy as np
import subprocess
import sys

from . import check
from . import context
from .obs import expect, simulate
from . import resample
from . import state
from . import version


def dtype_value(value, name='value'):
    """The dtype for columns that store observation values."""
    if hasattr(value, 'dtype'):
        value_dtype = value.dtype
    elif isinstance(value, int):
        value_dtype = np.dtype(int)
    elif isinstance(value, float):
        value_dtype = np.dtype(float)
    else:
        raise ValueError("no dtype for observed value {}".format(value))

    return (name, value_dtype)


def obs_types(params, obs_list, obs_reqd=False):
    """
    Return a sorted list of ``(unit, period)`` tuples that define the unique
    combinations of observation unit and observation period in the provided
    list of observations. If this list is empty, the observation models will
    instead be queried.

    :param params: The simulation parameters, or a simulation context.
    :type params: Union[dict, pypfilt.Context]
    :param obs_list: The list of observations for the simulation.
    :param obs_reqd: Whether a non-empty list of observations is required; if
        ``True``, an empty list will be return when no observations are
        provided.
    """
    if len(obs_list) > 0:
        return sorted(set((o['unit'], o['period']) for o in obs_list))
    elif obs_reqd:
        # Observations are required, return an empty list.
        return []
    elif isinstance(params, context.Context):
        # Query each observation model to determine the observation types.
        return sorted({(u, params.component['obs'][u].period)
                       for u in params.params['obs']})
    else:
        # Query each observation model to determine the observation types.
        return sorted({(u, params['component']['obs'][u].period)
                       for u in params['obs']})


class Table(abc.ABC):
    """
    The base class for summary statistic tables.

    Tables are used to record rows of summary statistics as a simulation
    progresses.
    """

    @abc.abstractmethod
    def dtype(self, ctx, obs_list, name):
        """
        Return the column names and data types, represented as a list of
        ``(name, data type)`` tuples. See the NumPy documentation for details.

        :param params: The simulation parameters.
        :param obs_list: A list of all observations.
        :param name: The table's name.
        """
        pass

    @abc.abstractmethod
    def n_rows(self, start_date, end_date, n_days, n_sys, forecasting):
        """
        Return the number of rows required for a single simulation.

        :param start_date: The date at which the simulation starts.
        :param end_date: The date at which the simulation ends.
        :param n_days: The number of days for which the simulation runs.
        :param n_sys: The number of observation systems (i.e., data sources).
        :param forecasting: ``True`` if this is a forecasting simulation,
            otherwise ``False``.
        """
        pass

    @abc.abstractmethod
    def add_rows(self, hist, weights, fs_date, dates, obs_types, insert_fn):
        """
        Record rows of summary statistics for some portion of a simulation.

        :param hist: The particle history matrix.
        :param weights: The weight of each particle at each date in the
            simulation window; it has dimensions ``(d, p)`` for ``d`` days and
            ``p`` particles.
        :param fs_date: The forecasting date; if this is not a forecasting
            simulation, this is the date at which the simulation ends.
        :param dates: A list of ``(datetime, ix, hist_ix)`` tuples that
            identify each day in the simulation window, the index of that day
            in the simulation window, and the index of that day in the
            particle history matrix.
        :param obs_types: A set of ``(unit, period)`` tuples that identify
            each observation system from which observations have been taken.
        :param insert_fn: A function that inserts one or more rows into the
            underlying data table; see the examples below.

        The row insertion function can be used as follows:

        .. code-block:: python

           # Insert a single row, represented as a tuple.
           insert_fn((x, y, z))
           # Insert multiple rows, represented as a list of tuples.
           insert_fn([(x0, y0, z0), (x1, y1, z1)], n=2)
        """
        pass

    def finished(self, hist, weights, fs_date, dates, obs_types, insert_fn):
        """
        Record rows of summary statistics at the end of a simulation.

        The parameters are as per :meth:`.add_rows`.

        Derived classes should only implement this method if rows must be
        recorded by this method; the provided method does nothing.
        """
        pass


class Monitor(abc.ABC):
    """
    The base class for simulation monitors.

    Monitors are used to calculate quantities that:

    * Are used by multiple Tables (i.e., avoiding repeated computation); or
    * Require a complete simulation for calculation (as distinct from Tables,
      which incrementally record rows as a simulation progresses).

    The quantities calculated by a Monitor can then be recorded by
    :meth:`.Table.add_rows` and/or :meth:`.Table.finished`.
    """

    def prepare(self, ctx, obs_list, name):
        """
        Perform any required preparation prior to a set of simulations.

        :param params: The simulation parameters.
        :param obs_list: A list of all observations.
        :param name: The monitor's name.
        """
        pass

    def begin_sim(self, start_date, end_date, n_days, n_sys, forecasting):
        """
        Perform any required preparation at the start of a simulation.

        :param start_date: The date at which the simulation starts.
        :param end_date: The date at which the simulation ends.
        :param n_days: The number of days for which the simulation runs.
        :param n_sys: The number of observation systems (i.e., data sources).
        :param forecasting: ``True`` if this is a forecasting simulation,
            otherwise ``False``.
        """
        pass

    @abc.abstractmethod
    def monitor(self, hist, weights, fs_date, dates, obs_types):
        """
        Monitor the simulation progress.

        :param hist: The particle history matrix.
        :param weights: The weight of each particle at each date in the
            simulation window; it has dimensions ``(d, p)`` for ``d`` days and
            ``p`` particles.
        :param fs_date: The forecasting date; if this is not a forecasting
            simulation, this is the date at which the simulation ends.
        :param dates: A list of ``(datetime, ix, hist_ix)`` tuples that
            identify each day in the simulation window, the index of that day
            in the simulation window, and the index of that day in the
            particle history matrix.
        :param obs_types: A set of ``(unit, period)`` tuples that identify
            each observation system from which observations have been taken.
        """
        pass

    def end_sim(self, hist, weights, fs_date, dates, obs_types):
        """
        Finalise the data as required for the relevant summary statistics.

        The parameters are as per :meth:`.monitor`.

        Derived classes should only implement this method if finalisation of
        the monitored data is required; the provided method does nothing.
        """
        pass

    @abc.abstractmethod
    def load_state(self, grp):
        """
        Load the monitor state from a cache file.

        :param grp: The h5py Group object from which to load the state.
        """
        pass

    @abc.abstractmethod
    def save_state(self, grp):
        """
        Save the monitor state to a cache file.

        :param grp: The h5py Group object in which to save the state.
        """
        pass


class ParamCovar(Table):
    """
    Calculate the covariance between all pairs of model parameters during each
    simulation.
    """

    def dtype(self, ctx, obs_list, name):
        self.__ctx = ctx
        self.__params = ctx.params
        # Only calculate covariances between model parameters that admit
        # continuous kernels.
        details = ctx.component['model'].describe()
        self.__param_info = [(n, ix)
                             for (ix, (n, smooth, _, _)) in enumerate(details)
                             if smooth]
        self.__param_names = [info[0] for info in self.__param_info]
        self.__num_params = len(self.__param_info)

        fs_date = ctx.component['time'].dtype('fs_date')
        date = ctx.component['time'].dtype('date')
        param1 = ('param1', h5py.string_dtype())
        param2 = ('param2', h5py.string_dtype())
        covar = ('covar', np.float64)
        return [fs_date, date, param1, param2, covar]

    def n_rows(self, start_date, end_date, n_days, n_sys, forecasting):
        num_params = len(self.__param_info)
        return n_days * num_params * (num_params - 1) // 2

    def add_rows(self, hist, weights, fs_date, dates, obs_types, insert_fn):
        from . import stats

        num_params = len(self.__param_info)
        fs_date_enc = self.__ctx.component['time'].to_dtype(fs_date)

        for date, ix, hist_ix in dates:
            date_enc = self.__ctx.component['time'].to_dtype(date)
            min_ix = min(info[1] for info in self.__param_info)
            max_ix = max(info[1] for info in self.__param_info)
            x = hist[hist_ix, :, min_ix:max_ix + 1]
            covars = stats.cov_wt(x, weights[ix, :])
            for ix1 in range(num_params):
                name1 = self.__param_names[ix1]
                for ix2 in range(ix1 + 1, num_params):
                    name2 = self.__param_names[ix2]
                    row = (fs_date_enc, date_enc, name1, name2,
                           covars[ix1, ix2])
                    insert_fn(row)


class ModelCIs(Table):
    """
    Calculate fixed-probability central credible intervals for all state
    variables and model parameters.

    :param probs: an array of probabilities that define the size of each
        central credible interval.
        The default value is ``numpy.uint8([0, 50, 90, 95, 99, 100])``.
    """

    def __init__(self, probs=None):
        if probs is None:
            probs = np.uint8([0, 50, 90, 95, 99, 100])
        self.__probs = probs

    def dtype(self, ctx, obs_list, name):
        self.__ctx = ctx
        self.__params = ctx.params
        details = ctx.component['model'].describe()
        self.__sv_info = [(info[0], ix) for (ix, info) in enumerate(details)]
        self.__stat_info = ctx.component['model'].stat_info()
        self.__num_stats = len(self.__sv_info) + len(self.__stat_info)
        self.__probs = np.uint8(ctx.params.get_chained(
            ['summary', 'tables', name, 'credible_intervals'],
            self.__probs))

        fs_date = ctx.component['time'].dtype('fs_date')
        date = ctx.component['time'].dtype('date')
        prob = ('prob', np.int8)
        ymin = ('ymin', np.float64)
        ymax = ('ymax', np.float64)
        # State variables/parameters ('model') or statistics ('stat').
        value_type = ('type', h5py.string_dtype())
        name = ('name', h5py.string_dtype())
        return [fs_date, date, prob, ymin, ymax, value_type, name]

    def n_rows(self, start_date, end_date, n_days, n_sys, forecasting):
        # Need a row for each interval, for each day, for each parameter,
        # variable and statistic.
        return len(self.__probs) * n_days * self.__num_stats

    def add_rows(self, hist, weights, fs_date, dates, obs_types, insert_fn):
        from . import stats

        fs_date_enc = self.__ctx.component['time'].to_dtype(fs_date)
        for date, ix, hist_ix in dates:
            date_enc = self.__ctx.component['time'].to_dtype(date)

            # Identify which state vectors to examine.
            valid = self.__ctx.component['model'].is_valid(hist[hist_ix])
            # Note that np.where() returns a *tuple*  of arrays (one for
            # each dimension) and we're only scanning a 1D array.
            mask = np.where(valid)[0]

            if np.count_nonzero(mask):
                ws = weights[ix, mask]
                for (val, vix) in self.__sv_info:
                    sub_hist = hist[hist_ix, mask, vix]
                    cred_ints = stats.cred_wt(sub_hist, ws, self.__probs)
                    for cix, pctl in enumerate(self.__probs):
                        row = (fs_date_enc, date_enc, pctl,
                               cred_ints[pctl][0], cred_ints[pctl][1],
                               'model', val)
                        insert_fn(row)
                for (val, stat_fn) in self.__stat_info:
                    stat_vec = stat_fn(hist[hist_ix, mask])
                    cred_ints = stats.cred_wt(stat_vec, ws, self.__probs)
                    for cix, pctl in enumerate(self.__probs):
                        row = (fs_date_enc, date_enc, pctl,
                               cred_ints[pctl][0], cred_ints[pctl][1],
                               'stat', val)
                        insert_fn(row)
            else:
                for pctl in self.__probs:
                    for (val, _) in self.__sv_info:
                        row = (fs_date_enc, date_enc, pctl, 0, 0,
                               'model', val)
                        insert_fn(row)
                    for (val, _) in self.__stat_info:
                        row = (fs_date_enc, date_enc, pctl, 0, 0, 'stat', val)
                        insert_fn(row)


class Obs(Table):
    """
    Record the basic details of each observation; the columns are: ``'unit',
    'period', 'source', 'date', 'value', 'incomplete', 'upper_bound'``.
    """

    def dtype(self, ctx, obs_list, name):
        self.__written = False
        self.__obs_list = obs_list
        self.__ctx = ctx
        self.__params = ctx.params
        obs_val = self.__obs_list[0]['value']
        dtype = [('unit', h5py.string_dtype()),
                 ('source', h5py.string_dtype()),
                 ('period', np.int8),
                 ctx.component['time'].dtype('date'),
                 dtype_value(obs_val),
                 ('incomplete', np.bool),
                 dtype_value(obs_val, name='upper_bound')]
        return dtype

    def n_rows(self, start_date, end_date, n_days, n_sys, forecasting):
        if self.__written:
            return 0
        else:
            return len(self.__obs_list)

    def add_rows(self, hist, weights, fs_date, dates, obs_types, insert_fn):
        # Add all of the rows in finished(), below.
        pass

    def finished(self, hist, weights, fs_date, dates, obs_types, insert_fn):
        for o in self.__obs_list:
            date = self.__ctx.component['time'].to_dtype(o['date'])
            inc = o.get('incomplete', False)
            upr = o.get('upper_bound', 0)
            row = (o['unit'], o['source'], o['period'], date,
                   o['value'], inc, upr)
            insert_fn(row)
        self.__written = True


class SimulatedObs(Table):
    """
    Record simulated observations for each particle.
    """

    def dtype(self, ctx, obs_list, name):
        self.__ctx = ctx
        self.__params = ctx.params
        self.__rnd = np.random.default_rng(ctx.params.get('prng_seed'))
        self.__sim = np.random.default_rng(ctx.params.get('prng_seed'))
        unit = ('unit', h5py.string_dtype(encoding='utf-8'))
        period = ('period', np.int8)
        fs_date = ctx.component['time'].dtype('fs_date')
        date = ctx.component['time'].dtype('date')
        value = ('value', np.float64)
        return [unit, period, fs_date, date, value]

    def n_rows(self, start_date, end_date, n_days, n_sys, forecasting):
        # NOTE: when forecasting, we only want to resample the particles once,
        # at the start of the forecast. This ensures that the simulated
        # observations reflect individual model trajectories. Note that this
        # is not possible during the estimation pass, because any observation
        # can trigger resampling.
        self.__sample_ixs = None
        self.__forecasting = forecasting
        # Need one row for each particle, for each day, for each data source.
        n_px = self.__params['hist']['px_count']
        return n_px * n_days * n_sys

    def add_rows(self, hist, weights, fs_date, dates, obs_types, insert_fn):
        fs_date_enc = self.__ctx.component['time'].to_dtype(fs_date)

        for date, ix, hist_ix in dates:
            for (unit, period) in obs_types:
                date_enc = self.__ctx.component['time'].to_dtype(date)
                curr = hist[hist_ix]
                n_back = self.__params['steps_per_unit'] * period
                prev = state.earlier_states(hist, hist_ix, n_back)
                exp_values = expect(self.__ctx, date, unit, period,
                                    prev, curr)

                if self.__params['hist']['px_count'] > 1:
                    # NOTE: resample the particles so that weights are uniform.
                    if self.__sample_ixs is None:
                        (sample_ixs, _weight) = resample.resample_weights(
                            weights[ix, :], self.__rnd)
                        if self.__forecasting:
                            self.__sample_ixs = sample_ixs
                    else:
                        sample_ixs = self.__sample_ixs
                else:
                    # NOTE: with only a single particle, we need to convert
                    # sim_values from a scalar to an array.
                    exp_values = np.array([exp_values])
                    sample_ixs = [0]

                sim_values = simulate(self.__ctx, date, unit, period,
                                      exp_values[sample_ixs], rng=self.__sim)
                if len(sample_ixs) == 1:
                    sim_values = [sim_values]
                for value in sim_values:
                    insert_fn((unit, period, fs_date_enc, date_enc, value))


class PredictiveCIs(Table):
    """
    Record fixed-probability central credible intervals for the observations.

    :param exp_obs_monitor: a :class:`pypfilt.summary.ExpectedObsMonitor`.
    :param probs: an array of probabilities that define the size of each
        central credible interval.
        The default value is ``numpy.uint8([0, 50, 95])``.
    """

    def __init__(self, exp_obs_monitor, probs=None):
        if probs is None:
            probs = np.uint8([0, 50, 95])
        self.__probs = np.sort(probs)
        self.__monitor_name = exp_obs_monitor
        qtl_prs = self.__probs.astype(float) / 100
        qtl_lwr = 0.5 - 0.5 * qtl_prs
        qtl_upr = 0.5 + 0.5 * qtl_prs
        self.__qtls = np.sort(np.unique(np.r_[qtl_lwr, qtl_upr]))

    def dtype(self, ctx, obs_list, name):
        self.__monitor = ctx.component['summary_monitor'][self.__monitor_name]
        self.__ctx = ctx
        self.__params = ctx.params
        unit = ('unit', h5py.string_dtype())
        period = ('period', np.int8)
        fs_date = ctx.component['time'].dtype('fs_date')
        date = ctx.component['time'].dtype('date')
        prob = ('prob', np.int8)
        ymin = ('ymin', np.float64)
        ymax = ('ymax', np.float64)
        return [unit, period, fs_date, date, prob, ymin, ymax]

    def n_rows(self, start_date, end_date, n_days, n_sys, forecasting):
        # Need a row for each interval, for each day, for each data source.
        return len(self.__probs) * n_days * n_sys

    def add_rows(self, hist, weights, fs_date, dates, obs_types, insert_fn):
        exp_obs = self.__monitor.expected_obs
        fs_date_enc = self.__ctx.component['time'].to_dtype(fs_date)

        for (unit, period) in obs_types:
            exp_sys = exp_obs[(unit, period)]
            op = self.__params['obs'][unit]
            obs_model = self.__ctx.component['obs'][unit]

            for date, ix, _ in dates:
                date_enc = self.__ctx.component['time'].to_dtype(date)
                mu = exp_sys[ix, :]
                ws = weights[ix, :]
                # Quantiles are ordered from smallest to largest.
                qtls = obs_model.quantiles(self.__ctx, op, date, mu, ws,
                                           self.__qtls)
                # Iterate from broadest to narrowest CI.
                for ix, pr in enumerate(self.__probs[::-1]):
                    row = (unit, period, fs_date_enc, date_enc, pr,
                           qtls[ix], qtls[- (ix + 1)])
                    insert_fn(row)


class ExpectedObsMonitor(Monitor):
    """
    Record expected observations for each particle.

    This is typically an expensive operation, and this monitor allows multiple
    summary tables to obtain these values without recalculating them.
    """

    expected_obs = None
    """
    The expected observation for each particle for the duration of the
    **current simulation window**.

    Note that this is **only** valid for tables to inspect in each call to
    ``add_rows()``, and **not** in a call to ``finished()``.
    """

    def __init__(self):
        self.expected_obs = None

    def prepare(self, ctx, obs_list, name):
        self.__ctx = ctx
        self.__params = ctx.params
        self.expected_obs = None

    def monitor(self, hist, weights, fs_date, dates, obs_types):
        """Record the peak for each particle during a forecasting run."""
        self.expected_obs = {}

        # Do nothing more if there are no dates to summarise.
        num_dates = len(dates)
        if num_dates == 0:
            # Ensure an empty data structure exists, at least.
            for (u, p) in obs_types:
                self.expected_obs[u, p] = np.array([])
            return

        periods = set([p for (_, p) in obs_types])

        times = [date for (date, ix, hist_ix) in dates]
        exp_shape = (len(times), self.__params['size'])
        for (u, p) in obs_types:
            self.expected_obs[u, p] = np.zeros(exp_shape)

        date_ix = 0
        for date, ix, hist_ix in dates:
            curr = hist[hist_ix]
            # Record the expected observations.
            for p in periods:
                n_back = self.__params['steps_per_unit'] * p
                prev = state.earlier_states(hist, hist_ix, n_back)

                valid_types = [(u, pd) for (u, pd) in obs_types if p == pd]
                for (u, p) in valid_types:
                    values = expect(self.__ctx, date, u, p, prev, curr)
                    self.expected_obs[u, p][date_ix] = values

            date_ix += 1

    def end_sim(self, hist, weights, fs_date, dates, obs_types):
        self.expected_obs = None

    def load_state(self, grp):
        """Load the monitor state for disk."""
        # NOTE: nothing needs to be loaded.
        pass

    def save_state(self, grp):
        """Save the monitor state to disk."""
        # NOTE: nothing needs to be saved.
        pass


class HDF5(object):
    """
    Save tables of summary statistics to an HDF5 file.

    :param params: The simulation parameters.
    :param obs_list: A list of all observations.
    """

    def __init__(self, params, obs_list):
        # Store simulation metadata.
        meta = Metadata()
        self.__metadata = meta.build(params)

        # Store the observations.
        self.__all_obs = obs_list

        # Allocate variables to store the details of each summary table.
        self.__tbl_dict = {}
        self.__dtypes = {}
        # When a simulation commences, this will be a dictionary that maps
        # table names to NumPy structured arrays; the value of ``None``
        # indicates that no tables have been allocated.
        self.__df = None

        self.__monitors = {}

        self.__data_group = 'data'

        # If True when self.__only_fs is True, the current simulation is not a
        # forecasting simulation, and tables should be ignored.
        # Note that monitors are *never* ignored.
        self.__ignore = False

        # Identify all combinations for observation periods and units.
        self.__obs_types = obs_types(params, obs_list)

    def initialise(self, ctx):
        """
        Initialise each table and monitor.

        This should be called before running a single set of simulations, and
        is called by each of the top-level pypfilt functions.
        """
        self.__ctx = ctx
        self.__params = ctx.params

        # If True, start statistics at the beginning of the simulation period,
        # otherwise start statistics at the date of the first observation.
        self.__first_day = ctx.params['summary']['from_first_day']
        # If True, only calculate statistics for forecasting simulations.
        self.__only_fs = ctx.params['summary']['only_forecasts']

        self.__monitors = {}
        for (name, monitor) in ctx.component['summary_monitor'].items():
            if name in self.__monitors:
                raise ValueError("Monitor '{}' already exists".format(name))
            self.__monitors[name] = monitor
            # NOTE: provide the monitor name here so that the monitor can
            # look for monitor-specific parameters.
            monitor.prepare(self.__ctx, self.__all_obs, name)

        self.__tbl_dict = {}
        for (name, table) in ctx.component['summary_table'].items():
            if name in self.__tbl_dict:
                raise ValueError("Table '{}' already exists".format(name))
            self.__tbl_dict[name] = table
            # NOTE: provide the table name here so that the table can look for
            # table-specific parameters.
            self.__dtypes[name] = table.dtype(self.__ctx, self.__all_obs,
                                              name)

    def load_state(self, grp):
        """
        Load the internal state of each monitor from a cache file.

        :param grp: The h5py Group object from which to load the state.
        """
        grp_names = set()
        for (name, mon) in self.__monitors.items():
            if name in grp_names:
                msg = "Multiple {} monitors".format(name)
                raise ValueError(msg)
            else:
                grp_names.add(name)
            mon_grp = grp.require_group(name)
            mon.load_state(mon_grp)

    def save_state(self, grp):
        """
        Save the internal state of each monitor to a cache file.

        :param grp: The h5py Group object in which to save the state.
        """
        grp_names = set()
        for (name, mon) in self.__monitors.items():
            if name in grp_names:
                msg = "Multiple {} monitors".format(name)
                raise ValueError(msg)
            else:
                grp_names.add(name)
            mon_grp = grp.require_group(name)
            mon.save_state(mon_grp)

    def allocate(self, ctx, start_date, end_date, forecasting=False):
        """
        Allocate space for the simulation statistics.

        This is called by ``pypfilt.pfilter.run`` before running a simulation.
        When multiple simulations are run, such as a series of forecasts, this
        will be called before each simulation.
        """
        if self.__df is not None:
            raise ValueError("Tables have already been allocated")

        if self.__only_fs and not forecasting:
            # Flag this simulation as being ignored.
            self.__ignore = True
        else:
            self.__ignore = False

        logger = logging.getLogger(__name__)

        self.__start_date = start_date
        self.__end_date = end_date
        if forecasting:
            # Forecasting from the start of the simulation period.
            self.__fs_date = start_date
        else:
            # Not forecasting, so all observations are included.
            # NOTE: set this to the end of the *true* simulation period.
            # This avoids a problem where, when using a cache file, the
            # estimation run will end at the final forecasting data and both
            # the estimation run and this final forecast will be identified by
            # the same "forecasting" date.
            self.__fs_date = ctx.params['time']['until']

        # Retain only observations that occur during the simulation period.
        def include_obs(d):
            try:
                return start_date <= d <= end_date
            except TypeError:
                raise ValueError("invalid observation date '{}'".format(d))

        fs_obs = [o for o in self.__all_obs if include_obs(o['date'])]

        # Calculate the number of days from the first observation date to the
        # end of the simulation
        obs_dates = sorted(set(o['date'] for o in fs_obs))
        # Either use observation times *or* each time unit.
        steps_per_unit = self.__params['steps_per_unit']
        if self.__first_day:
            t0 = start_date
        else:
            t0 = min(obs_dates)
        self.__all_dates = [t0]
        for (step, time) in self.__ctx.component['time'].steps():
            if step % steps_per_unit == 0 and time > t0:
                self.__all_dates.append(time)
        n_days = len(self.__all_dates)
        n_sys = len(self.__obs_types)

        # Notify each monitor, regardless of whether tables are ignored.
        for mon in self.__monitors.values():
            mon.begin_sim(start_date, end_date, n_days, n_sys, forecasting)

        if self.__ignore:
            logger.debug("Summary.allocate({}, {}): {} points".format(
                start_date, end_date, n_days))
            return

        # Determine the number of rows to allocate for each table.
        # Note: the items() method returns a list in Python 2 and a view
        # object in Python 3; since the number of tables will always be small
        # (on the order of 10) the overhead of using items() in Python 2 and
        # not iteritems() --- which returns an interator --- is negligible.
        n_rows = {n: t.n_rows(start_date, end_date, n_days, n_sys,
                              forecasting)
                  for (n, t) in self.__tbl_dict.items()}
        # Allocate tables that require at least one row.
        self.__df = {tbl: np.empty(n_rows[tbl], dtype=self.__dtypes[tbl])
                     for tbl in n_rows if n_rows[tbl] > 0}
        # Initialise a row counter for each allocated table.
        self.__ix = {tbl: 0 for tbl in self.__df}
        # Create a row insertion function for each allocated table.
        self.__insert = {tbl: self.__insert_row_fn(tbl) for tbl in self.__df}

        logger.debug("Summary.allocate({}, {}): {} points".format(
            start_date, end_date, n_days))

    def __insert_row_fn(self, tbl):
        def insert(fields, n=1):
            row_ix = self.__ix[tbl]
            self.__df[tbl][row_ix:row_ix + n] = fields
            self.__ix[tbl] += n
        return insert

    def summarise(self, hist, sim_time, start_date, end_date, offset):
        """
        Calculate statistics for some portion of the simulation period.

        :param hist: The particle history matrix.
        :param sim_time: The entire simulation period.
        :param start_date: The start of the period to be summarised.
        :param end_date: The end of the period to be summarised.
        :param offset: The history matrix time-step offset.
        """
        if self.__df is None and not self.__ignore:
            raise ValueError("Tables have not been allocated")

        if self.__start_date > start_date or self.__end_date < end_date:
            raise ValueError("Summary.summarise() called for invalid period")

        # Ensure we have received the entire particle history matrix.
        check.is_entire_matrix(self.__params, hist)

        logger = logging.getLogger(__name__)

        dates = [(d, ix, sim_time.step_of(d) + offset)
                 for (ix, d) in enumerate(d for d in self.__all_dates
                                          if start_date <= d <= end_date)]
        num_dates = len(dates)

        logger.debug("Summary.summarise({}, {}): {} dates".format(
            start_date, end_date, num_dates))

        # Record the particle weights.
        weights = np.zeros((num_dates, hist.shape[1]))
        for date, ix, hist_ix in dates:
            weights[ix, :] = hist[hist_ix, :, -2]

        fs_date = self.__fs_date
        obs_types = self.__obs_types

        for mon in self.__monitors.values():
            mon.monitor(hist, weights, fs_date, dates, obs_types)

        tables = self.__df if self.__df is not None else []

        for tbl in tables:
            insert_fn = self.__insert[tbl]
            self.__tbl_dict[tbl].add_rows(
                hist, weights, fs_date, dates, obs_types, insert_fn)

        if end_date == self.__end_date:
            for mon in self.__monitors.values():
                mon.end_sim(hist, weights, fs_date, dates, obs_types)

            for tbl in tables:
                insert_fn = self.__insert[tbl]
                self.__tbl_dict[tbl].finished(
                    hist, weights, fs_date, dates, obs_types, insert_fn)

    def get_stats(self):
        """Return the calculated statistics for a single simulation."""
        if self.__df is None:
            if self.__ignore:
                return {}
            else:
                raise ValueError("Tables have not been created")

        logger = logging.getLogger(__name__)
        logger.debug("Summary.get()")

        # Check all table rows are filled (and no further).
        for tbl in self.__df:
            alloc = self.__df[tbl].shape[0]
            used = self.__ix[tbl]
            if alloc != used:
                msg = "Table '{}' allocated {} rows but filled {}"
                raise ValueError(msg.format(tbl, alloc, used))

        # Return the summary tables and remove them from this class instance.
        stats = self.__df
        self.__df = None
        return stats

    def save_forecasts(self, fs, filename):
        """
        Save forecast summaries to disk in the HDF5 binary data format.

        This function creates the following datasets that summarise the
        estimation and forecasting outputs:

        - ``'data/TABLE'`` for each table.

        The provided metadata will be recorded under ``'meta/'``.

        If dataset creation timestamps are enabled, two simulations that
        produce identical outputs will not result in identical files.
        Timestamps will be disabled where possible (requires h5py >= 2.2):

        - ``'hdf5_track_times'``: Presence of creation timestamps.

        :param fs: Simulation outputs, as returned by ``pypfilt.forecast()``.
        :param filename: The filename to which the data will be written.
        """
        fs_dates = [d for d in fs.keys() if d not in ['obs', 'complete']]
        fs_compl = fs_dates[:]
        # Check that the estimation pass results exist.
        if 'complete' in fs.keys():
            fs_compl.append('complete')

        # Construct aggregate data tables.
        # Note that some tables may not exist for every simulation.
        tbl_dict = {}
        tbl_names = set(ns for fd in fs_compl for ns in fs[fd]['summary'])
        for n in tbl_names:
            in_fs = [fdate for fdate in fs_compl
                     if n in fs[fdate]['summary']]
            tbl_dict[n] = np.concatenate([fs[fdate]['summary'][n]
                                          for fdate in in_fs])

        # Attempt to avoid tracking times (not supported by h5py < 2.2).
        kwargs = {'track_times': False}

        def save_data(g, name, value):
            if isinstance(value, dict):
                sub_g = g.create_group(name)
                for sub_key, sub_value in sorted(value.items()):
                    save_data(sub_g, sub_key, sub_value)
            else:
                try:
                    g.create_dataset(name, data=value, **kwargs)
                except TypeError:
                    msg = 'Error saving dataset "{}" with value {} and type {}'
                    raise ValueError(msg.format(name, value,
                                                type(value).__name__))

        with h5py.File(filename, 'w') as f:
            # Record whether tracking times have been disabled.
            try:
                f.create_dataset('hdf5_track_times', data=False, **kwargs)
            except TypeError:
                # Cannot disable tracking times (h5py < 2.2).
                kwargs = {}
                f.create_dataset('hdf5_track_times', data=True, **kwargs)

            # Save the associated metadata, if any.
            if self.__metadata:
                meta_grp = f.create_group('meta')
                for k, v in sorted(self.__metadata.items()):
                    save_data(meta_grp, k, v)

            # Compress and checksum the data tables.
            kwargs['compression'] = 'gzip'
            kwargs['shuffle'] = True
            kwargs['fletcher32'] = True

            # Save the data tables.
            df_grp = f.create_group(self.__data_group)
            for tbl in tbl_dict:
                df_grp.create_dataset(tbl, data=tbl_dict[tbl], **kwargs)


class FakeVal(object):
    """A fake value that converts basic arithmetic operations into strings."""

    def __init__(self, value):
        self.value = value

    def __rtruediv__(self, x):
        return FakeVal("({} / {})".format(x, self.value))

    def __str__(self):
        return self.value

    def __repr__(self):
        return repr(self.value)


class FakeRNG(object):
    """
    A fake random number generator whose methods return strings that describe
    each method call.
    """

    def __init__(self, prefix="random."):
        self.__prefix = prefix

    def __getattr__(self, name):
        """Called when an attribute lookup has not found the attribute."""
        name = self.__prefix + name

        def log_call(*args, **kwargs):
            """Return a string representation of a method call."""
            # Ignore the 'size' keyword argument if its value is None.
            kw_iter = [(k, v) for (k, v) in kwargs.items()
                       if k != 'size' or v is not None]
            return FakeVal("{}({})".format(
                name,
                ", ".join(
                    [repr(a) for a in args] +
                    ["{}={}".format(k, repr(v)) for k, v in kw_iter]
                )))

        return log_call


class Metadata(object):
    """
    Document the simulation parameters and system environment for a set of
    simulations. A black-list (``ignore_dict``) defines which members of
    the parameters dictionary will be excluded from this metadata, see
    :func:`~Metadata.filter` for details.
    """

    def __init__(self):
        # Ignore the following parameters; some will be replaced by better
        # descriptions (e.g., 'prior') and some should always be ignored
        # because they have no meaningful representation.
        self.ignore_dict = {
            'model': {
                'prior': None,
            },
            'hist': {
                'extra_col_fns': None,
            },
            'component': {
                None: None,
            },
            'data': None,
            'hooks': None,
            'px_range': None,
        }

    def build(self, params):
        """
        Construct a metadata dictionary that documents the simulation
        parameters and system environment. Note that this should be generated
        at the **start** of the simulation, and that the git metadata will
        only be valid if the working directory is located within a git
        repository.

        :param params: The simulation parameters.

        By default, the versions of ``pypfilt``, ``h5py``, ``numpy`` and
        ``scipy`` are recorded. The following example demonstrates how to also
        record the installed version of the ``epifx`` package:

        .. code-block:: python

           import epifx
           import pypfilt.summary
           params = ...
           params['summary']['meta']['packages'].append('epifx')
           meta = pypfilt.summary.Metadata()
        """

        # Import modules for extracting system-specific details.
        import platform
        # Import modules for extracting package versions.
        import scipy
        # Import modules for parsing scenario configurations.
        import toml

        # Record the command line used to launch this simulation.
        # Note that sys.argv is a list of native strings.
        cmdline = " ".join(sys.argv)

        meta = {
            'package': {
                'python': platform.python_version(),
                'h5py': self.pkg_version(h5py),
                'pypfilt': self.pkg_version(version),
                'numpy': self.pkg_version(np),
                'scipy': self.pkg_version(scipy),
                'toml': self.pkg_version(toml),
            },
            'sim': {
                'cmdline': cmdline,
            },
            'param': self.filter(params, self.ignore_dict, self.encode),
            'prior': self.priors(params),
        }

        git_data = self.git_data()
        if git_data:
            meta['git'] = git_data

        # Record the versions of user-specified packages (if any).
        pkgs = params.get('summary', {}).get('meta', {}).get('packages', [])
        for package_name in pkgs:
            mod_obj = importlib.import_module(package_name)
            meta['package'][package_name] = self.pkg_version(mod_obj)

        # Record all of the component object names.
        meta['param']['component'] = self.object_names(params['component'])

        return meta

    def filter(self, values, ignore, encode_fn):
        """
        Recursively filter items from a dictionary, used to remove parameters
        from the metadata dictionary that, e.g., have no meaningful
        representation.

        :param values: The original dictionary.
        :param ignore: A dictionary that specifies which values to ignore.
        :param encode_fn: A function that encodes the remaining values (see
            :func:`.encode_value`).

        For example, to ignore ``['px_range']``,  ``['random']``, and
        ``'obs_model'`` and ``'obs_llhd'`` for *every* observation system
        when using ``epifx``:

        .. code-block:: python

           m = pypfilt.summary.Metadata()
           ignore = {
               'px_range': None,
               'random': None,
               # Note the use of ``None`` to match any key under 'obs'.
               'obs': {None: {'obs_model': None, 'obs_llhd': None}}
           }
           m.filter(params, ignore, m.encode)
        """
        retval = {}
        for k, v in values.items():
            if k in ignore and ignore[k] is None:
                # Ignore this parameter.
                continue
            elif k in ignore and isinstance(v, dict):
                # Recursively encode this dictionary, maybe ignoring values.
                retval[k] = self.filter(v, ignore[k], encode_fn)
            elif None in ignore and isinstance(v, dict):
                if ignore[None] is None:
                    continue
                # Recursively encode this dictionary, maybe ignoring values.
                retval[k] = self.filter(v, ignore[None], encode_fn)
            elif isinstance(v, dict):
                # Recursively encode this dictionary without ignoring values.
                retval[k] = self.filter(v, {}, encode_fn)
            else:
                retval[k] = encode_fn(v)
        return retval

    def encode(self, value):
        """
        Encode values in a form suitable for serialisation in HDF5 files.

        * Integer values are converted to ``numpy.int32`` values.
        * Floating-point values and arrays retain their data type.
        * All other (i.e., non-numerical) values are converted to UTF-8
          strings.
        """
        if isinstance(value, (int, np.int64)):
            # Avoid storing 64-bit integers since R doesn't support them.
            return np.int32(value)
        elif isinstance(value, (float, np.ndarray)):
            # Ensure that numerical values retain their data type.
            return value
        elif isinstance(value, str):
            return value
        else:
            # Convert non-numerical values to UTF-8 strings.
            return str(value)

    def object_names(self, object_dict):
        """
        Return the fully qualified name of each object in a (possibly nested)
        dictionary.
        """
        names = {}
        for (name, value) in object_dict.items():
            if isinstance(value, dict):
                value_dict = self.object_names(value)
                if value_dict:
                    names[name] = value_dict
            else:
                names[name] = self.object_name(value)
        return names

    def object_name(self, obj):
        """
        Return the fully qualified name of the object as a byte string.
        """
        if hasattr(obj, '__name__'):
            # The object is a class.
            obj_name = obj.__name__
            obj_mod = obj.__module__
        elif hasattr(obj, '__class__'):
            # The object is a class instance.
            obj_name = obj.__class__.__name__
            obj_mod = obj.__class__.__module__
        else:
            # The object doesn't appear to be a class or a class instance.
            obj_name = str(obj)
            obj_mod = None

        if obj_mod is not None:
            descr = "{}.{}".format(obj_mod, obj_name)
        else:
            descr = obj_name

        return descr

    def priors(self, params):
        """
        Return a dictionary that describes the model parameter priors.

        Each key identifies a parameter (by name); the corresponding value is
        a byte string representation of the prior distribution, which is
        typically a ``numpy.random.Generator`` method call.

        For example:

        .. code-block:: python

           {'R0': b'random.uniform(1.0, 2.0)',
            'gamma': b'(1 / random.uniform(1.0, 3.0))'}
        """
        logger = logging.getLogger(__name__)
        priors = params['model']['prior'].items()
        rng = FakeRNG()
        descr = {}
        for k, v in priors:
            try:
                descr[k] = str(v(rng))
            except TypeError:
                # Can occur when performing arbitrary operations.
                logger.debug("Can not describe the prior for '{}'".format(k))
                descr[k] = "unknown"
        return descr

    def pkg_version(self, module):
        """Attempt to obtain the version of a Python module."""
        try:
            return str(module.__version__)
        except AttributeError:
            try:
                # Older versions of h5py store the version number here.
                return str(module.version.version)
            except AttributeError:
                return "unknown"

    def git_data(self):
        """
        Record the status of the git repository within which the working
        directory is located (if such a repository exists).
        """
        # Determine the encoding for the default locale.
        default_encoding = locale.getdefaultlocale()[1]
        logger = logging.getLogger(__name__)
        enc_msg = "Extracting git metadata using locale encoding '{}'"
        logger.debug(enc_msg.format(default_encoding))

        # Return no metadata if git is not installed, or if the working
        # directory is not located within a git repository.
        try:
            git_head = self.run_cmd(['git', 'rev-parse', 'HEAD'])
        except FileNotFoundError:
            logger.info('Could not run git; do you have git installed?')
            return {}
        if not git_head:
            logger.info('No HEAD commit; presumably not in a git repository?')
            return {}

        git_branch = self.run_cmd(['git', 'symbolic-ref', '--short', 'HEAD'])
        git_mod_files = self.run_cmd(['git', 'ls-files', '--modified'],
                                     all_lines=True, err_val=[])
        git_mod_files.sort()
        git_mod = len(git_mod_files) > 0
        return {
            'HEAD': git_head,
            'branch': git_branch,
            'modified': git_mod,
            'modified_files': [f.encode("utf-8") for f in git_mod_files],
        }

    def run_cmd(self, args, all_lines=False, err_val=''):
        """
        Run a command and return the (Unicode) output. By default, only the
        first line is returned; set ``all_lines=True`` to receive all of the
        output as a list of Unicode strings. If the command returns a non-zero
        exit status, return ``err_val`` instead.
        """
        try:
            # Return the output as a single byte string.
            lines = subprocess.check_output(args, stderr=subprocess.STDOUT)
        except subprocess.CalledProcessError:
            return err_val
        # Decode and break into lines according to Unicode boundaries.
        # For details see:
        # * https://docs.python.org/2/library/stdtypes.html#unicode.splitlines
        # * https://docs.python.org/3/library/stdtypes.html#str.splitlines
        default_encoding = locale.getdefaultlocale()[1]
        lines = lines.decode(default_encoding).splitlines()
        if all_lines:
            return lines
        else:
            return lines[0]


def obs_table(ctx, all_obs):
    """
    Create a structured array that contains the details of each observation.
    The columns are: ``'unit', 'period', 'source', 'date', 'value',
    'incomplete', 'upper_bound'``.

    :param params: The simulation parameters.
    :param all_obs: A list of all observations.
    """
    unit = ('unit', h5py.string_dtype())
    period = ('period', np.int8)
    source = ('source', h5py.string_dtype())
    date = ctx.component['time'].dtype('date')
    obs_val = all_obs[0]['value']
    value = dtype_value(obs_val)
    incomplete = ('incomplete', np.bool)
    upper = dtype_value(obs_val, name='upper_bound')
    dtype = [unit, period, source, date, value, incomplete, upper]
    obs_df = np.empty(len(all_obs), dtype=dtype)
    for ix, o in enumerate(all_obs):
        date_enc = ctx.component['time'].to_dtype(o['date'])
        inc = o.get('incomplete', False)
        upr = o.get('upper_bound', 0)
        obs_df[ix] = (o['unit'], o['period'], o['source'],
                      date_enc, o['value'], inc, upr)
    return obs_df


def convert_cols(data, converters):
    """
    Convert columns in a structured array from one type to another.

    :param data: The input structured array.
    :param converters: A dictionary that maps (unicode) column names to
        ``(convert_fn, new_dtype)`` tuples, which contain a conversion
        function and define the output dtype.
    :returns: A new structured array.
    """
    dt = data.dtype
    col_names = dt.names
    conv_names = [n for n in converters.keys()]
    conv_fns = {}
    new_dtype = []
    for col in col_names:
        if col in conv_names:
            convert_fn, convert_dt = converters[col]
            conv_fns[col] = convert_fn
            new_dtype.append((col, convert_dt))
        else:
            new_dtype.append((col, dt.fields[col][0]))
    # Create and fill the new array.
    data_out = np.zeros(data.shape, dtype=new_dtype)
    for col in col_names:
        if col in conv_names:
            data_out[col] = conv_fns[col](data[col])
        else:
            data_out[col] = data[col]
    return data_out


def default_converters(time_scale):
    """
    Return a dictionary for converting the ``'fs_date'`` and ``'date'``
    columns from (see :func:`~convert_cols`).
    """
    def convert_fn(np_col):
        return np.array([time_scale.from_dtype(val) for val in np_col])
    return {
        'fs_date': (convert_fn, time_scale.native_dtype()),
        'date': (convert_fn, time_scale.native_dtype()),
    }
