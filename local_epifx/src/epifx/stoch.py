"""Stochastic compartmental models."""

import logging
from pathlib import Path
import numpy as np
import pypfilt
from .model import Model


class SEEIIR(Model):
    """A stochastic SEEIIR compartment model."""

    __info = [
        ("S_U", False, 0, 1),
        ("E1_U", False, 0, 1), ("E2_U", False, 0, 1),
        ("I1_U", False, 0, 1), ("I2_U", False, 0, 1), 
        ("R_U", False, 0, 1),

        ("S_V", False, 0, 1),
        ("E1_V", False, 0, 1), ("E2_V", False, 0, 1),
        ("I1_V", False, 0, 1), ("I2_V", False, 0, 1), 
        ("R_V", False, 0, 1),

        ("R0", False, 1.0, 2.0),
        ("sigma", True, 1/3, 2.0),
        ("gamma", True, 1/3, 1.0),
        ("t0", False, 0, 100),
        ("R0_ix", False, 0, 1e6),
        ("R0_val", False, 0, 100),
        ("adjustment", False, 0, 1),
        ("mean_Ei", False, 0, 1), ("mean_Et", False, 0, 1)]

    ix_S_U = 0
    ix_E1_U = 1
    ix_E2_U = 2
    ix_I1_U = 3
    ix_I2_U = 4
    ix_R_U = 5

    ix_S_V = 6
    ix_E1_V = 7
    ix_E2_V = 8
    ix_I1_V = 9
    ix_I2_V = 10
    ix_R_V = 11


    ix_R0 = 12
    ix_sigma = 13
    ix_gamma = 14
    ix_t0 = 15
    ix_R0_ix = 16
    ix_R0_val = 17

    ix_adjustment = 18

    ix_mean_Ei = 19
    ix_mean_Et = 20

    R0_order_map = np.arange(0, 1000, 1)

    def __init__(self):
        self.__R0_lookup = None
        self.__external_lookup = None
        self.__regularise_R0_ix = False

    def state_size(self):
        """Return the size of the state vector."""
        return len(self.__info)

    def population_size(self):
        return self.popn_size

    def init(self, ctx, vec):
        """Initialise a state vector.

        :param ctx: The simulation context.
        :param vec: An uninitialised state vector of correct dimensions (see
            :py:func:`~state_size`).
        """
        self.popn_size = ctx.params['model']['population_size']
        self.__R0_lookup = None
        self.__external_lookup = None
        self.__vaccinations_lookup = None
        self.__regularise_R0_ix = ctx.params.get_chained(
            ['model', 'regularisation', 'R0_ix'], default=False)

        prior = ctx.params['model']['prior']
        rnd_size = vec[..., 0].shape
        rnd = ctx.component['random']['model']

        num_exps = 10.0
        vec[..., :] = 0

        vec[..., self.ix_S_U] = self.popn_size - num_exps
        vec[..., self.ix_I1_U] = num_exps
        
        vec[..., self.ix_R0] = prior['R0'](rnd, size=rnd_size)
        vec[..., self.ix_sigma] = prior['sigma'](rnd, size=rnd_size)
        vec[..., self.ix_gamma] = prior['gamma'](rnd, size=rnd_size)
        vec[..., self.ix_t0] = prior['t0'](rnd, size=rnd_size)

        vec[..., self.ix_adjustment] = 1

        self.load_samples_file(ctx, vec)
        self.load_lookup_tables(ctx)
        self.init_lookup_values(ctx, vec)

    def sample_columns(self):
        """Identify parameters that can be saved and loaded."""
        ix_tbl = {
            'R0': self.ix_R0,
            'sigma': self.ix_sigma,
            'gamma': self.ix_gamma,
            't0': self.ix_t0,
            'R0_ix': self.ix_R0_ix,
            'adjustment': self.ix_adjustment
        }
        return ix_tbl

    def load_samples_file(self, ctx, vec):
        """Load initial parameter values from an external data file."""
        if 'prior_samples' not in ctx.params['model']:
            return

        logger = logging.getLogger(__name__)
        samples = ctx.params['model']['prior_samples']
        data_dir = Path(ctx.params['data_dir'])
        data_file = data_dir / samples['file']
        columns = [(name, np.float) for name in samples['columns']]

        tbl = pypfilt.io.read_table(data_file, columns)
        if tbl.shape != vec[..., 0].shape:
            raise ValueError('Incompatible data shapes: {} and {}'.format(
                vec[..., 0].shape, tbl.shape))

        ix_tbl = self.sample_columns()
        for name in samples['columns']:
            if name not in ix_tbl:
                raise ValueError('Unknown parameter {}'.format(name))

            ix = ix_tbl[name]
            vec[..., ix] = tbl[name]

            # NOTE: warn if sampled values exceed the parameter bounds.
            min_val = np.min(tbl[name])
            max_val = np.max(tbl[name])
            if min_val < ctx.params['model']['param_min'][ix]:
                logger.warning('Sampled value for {} outside bounds'
                               .format(name))
            elif max_val > ctx.params['model']['param_max'][ix]:
                logger.warning('Sampled value for {} outside bounds'
                               .format(name))
            # Clip the sampled values to enforce the parameter bounds.
            # The alternative is to leave the sample values as provided, in
            # which case they will only be clipped if post-regularisation is
            # enabled and the particles are resampled.
            vec[..., ix] = np.clip(vec[..., ix],
                                   ctx.params['model']['param_min'][ix],
                                   ctx.params['model']['param_max'][ix])

    def resume_from_cache(self, ctx):
        """
        A simulation will begin from a saved state, so the model must check
        whether any lookup tables are defined.

        :param ctx: The simulation context.
        """
        self.load_lookup_tables(ctx)

    def load_lookup_tables(self, ctx):
        """
        Allow R0 and imported cases to be provided via lookup tables.

        :param ctx: The simulation context.
        """
        logger = logging.getLogger(__name__)
        tables = ctx.component.get('lookup', {})

        # Check for the R0 lookup table.
        if 'R0' in tables:
            self.__R0_lookup = tables['R0']
            logger.info('Using lookup table for R0 with {} values'.format(
                self.__R0_lookup.value_count()))

        # Check for the external exposures lookup table.
        exp_table = 'external_exposures'
        if exp_table in tables:
            self.__external_lookup = tables[exp_table]
            logger.info(
                'Using lookup table for external exposures with {} values'
                .format(self.__external_lookup.value_count()))

        vacc_table = 'vaccinations'
        if vacc_table in tables:
            self.__vaccinations_lookup = tables[vacc_table]
            logger.info(
                'Using lookup table for vaccinations with {} values'
                .format(self.__vaccinations_lookup.value_count()))

    def init_lookup_values(self, ctx, vec):
        """
        Initialise the ``R0_ix`` values if an R0 lookup table is defined.

        :param ctx: The simulation context.
        :param vec: An uninitialised state vector of correct dimensions (see
            :py:func:`~state_size`).
        """
        if self.__R0_lookup is not None:
            num_values = self.__R0_lookup.value_count()
            if num_values > 1:
                rnd = ctx.component['random']['model']
                rnd_size = vec[..., 0].shape
                vec[..., self.ix_R0_ix] = rnd.integers(num_values,
                                                       size=rnd_size)
            else:
                vec[..., self.ix_R0_ix] = 0


    def update(self, ctx, step_date, dt, is_fs, prev, curr):
        """Perform a single time-step.

        :param ctx: The simulation context.
        :param step_date: The date and time of the current time-step.
        :param dt: The time-step size (days).
        :param is_fs: Indicates whether this is a forecasting simulation.
        :param prev: The state before the time-step.
        :param curr: The state after the time-step (destructively updated).
        """

        rnd = ctx.component['random']['model']
        params = ctx.params

        # Update parameters and lookup tables that are defined in self.init()
        # and which will not exist if we are resuming from a cached state.
        self.popn_size = ctx.params['model']['population_size']

        # Extract each parameter.
        R0 = prev[..., self.ix_R0].copy()
        sigma = prev[..., self.ix_sigma].copy()
        gamma = prev[..., self.ix_gamma].copy()
        t0 = prev[..., self.ix_t0].copy()
        R0_ix = np.around(prev[..., self.ix_R0_ix]).astype(int)
        adjustment = prev[..., self.ix_adjustment].copy()

        if self.__R0_lookup is not None:
            R0 = self.get_R0_from_lookup(ctx, step_date, is_fs, R0_ix, params)

        # Extract each compartment.
        S_U = prev[...,   self.ix_S_U].astype(int)
        E1_U = prev[..., self.ix_E1_U].astype(int)
        E2_U = prev[..., self.ix_E2_U].astype(int)
        I1_U = prev[..., self.ix_I1_U].astype(int)
        I2_U = prev[..., self.ix_I2_U].astype(int)
        R_U = prev[..., self.ix_R_U].astype(int)

        S_V = prev[...,   self.ix_S_V].astype(int)
        E1_V = prev[..., self.ix_E1_V].astype(int)
        E2_V = prev[..., self.ix_E2_V].astype(int)
        I1_V = prev[..., self.ix_I1_V].astype(int)
        I2_V = prev[..., self.ix_I2_V].astype(int)
        R_V = prev[..., self.ix_R_V].astype(int)



        external = self.get_external_exposure_from_lookup(step_date, R0.shape)
        vacc_rate, mean_vacc_Ei, mean_vacc_Et = self.get_vaccinations_from_lookup(step_date, R0.shape)

        tau_V = 1 - mean_vacc_Et
        xi_V = 1 - mean_vacc_Ei

        val_lambda = (I1_U + I2_U) + tau_V * (I1_V + I2_V)

        # Only calculate adj. factor for backcasts (expecting Reff to be held constant across forecasting period)
        if not is_fs:
            n = self.popn_size
            
            # Calculating our adjustment factor
            # This needs to be bettered verified in context with vaccination
            # Argument in terms on dividing by Iu + Iv \tau_V??
            denom = (S_U + S_V).astype(float)

            adjustment = np.divide(n, denom, out = np.zeros_like(denom), where = denom != 0)

            # Is this really necessary? Not sure!
            adjustment = np.nan_to_num(adjustment, nan = 0, posinf = 0)


        epoch = ctx.component['time'].to_scalar(ctx.params['epoch'])
        curr_t = ctx.component['time'].to_scalar(step_date)
        zero_mask = t0 > (curr_t - epoch)
        R0[zero_mask] = 0
        sigma[zero_mask] = 0
        gamma[zero_mask] = 0

        beta = R0 * adjustment * gamma

        # Calculate the rates at which an individual leaves each compartment.
        s_U_out_rate = dt * (beta * val_lambda + external) / self.popn_size
        s_V_out_rate = dt * (beta * val_lambda * xi_V) / self.popn_size


        e_out_rate = dt * 2 * sigma
        i_out_rate = dt * 2 * gamma

        # Sample the outflow rate for each compartment.
        s_U_out = rnd.binomial(  S_U, - np.expm1(- s_U_out_rate))
        e1_U_out = rnd.binomial(E1_U, - np.expm1(- e_out_rate))
        e2_U_out = rnd.binomial(E2_U, - np.expm1(- e_out_rate))
        i1_U_out = rnd.binomial(I1_U, - np.expm1(- i_out_rate))
        i2_U_out = rnd.binomial(I2_U, - np.expm1(- i_out_rate))


        s_V_out = rnd.binomial(  S_V, - np.expm1(- s_V_out_rate))
        e1_V_out = rnd.binomial(E1_V, - np.expm1(- e_out_rate))
        e2_V_out = rnd.binomial(E2_V, - np.expm1(- e_out_rate))
        i1_V_out = rnd.binomial(I1_V, - np.expm1(- i_out_rate))
        i2_V_out = rnd.binomial(I2_V, - np.expm1(- i_out_rate))

        # Calculate vaccinations for each unvaccinated group:
        n_U = (S_U + E1_U + E2_U + I1_U + I2_U + R_U).astype(float)
        vacc_increase = np.divide(vacc_rate, n_U, out=np.zeros_like(n_U), where = n_U != 0) * dt

        # Calculate movement of vaccinated individuals after accounting for SEEIIR processes
        s_U_to_V  = rnd.binomial(S_U -  s_U_out,  -np.expm1(-vacc_increase))
        e1_U_to_V = rnd.binomial(E1_U - e1_U_out, -np.expm1(-vacc_increase))
        e2_U_to_V = rnd.binomial(E2_U - e2_U_out, -np.expm1(-vacc_increase))
        i1_U_to_V = rnd.binomial(I1_U - i1_U_out, -np.expm1(-vacc_increase))
        i2_U_to_V = rnd.binomial(I2_U - i2_U_out, -np.expm1(-vacc_increase))
        r_U_to_V  = rnd.binomial(R_U,             -np.expm1(-vacc_increase))


        # Update the compartment values.
        # Unvaccinateds:
        curr[...,  self.ix_S_U] =  S_U -  s_U_out            - s_U_to_V
        curr[..., self.ix_E1_U] = E1_U +  s_U_out - e1_U_out - e1_U_to_V
        curr[..., self.ix_E2_U] = E2_U + e1_U_out - e2_U_out - e2_U_to_V
        curr[..., self.ix_I1_U] = I1_U + e2_U_out - i1_U_out - i1_U_to_V
        curr[..., self.ix_I2_U] = I2_U + i1_U_out - i2_U_out - i2_U_to_V
        curr[..., self.ix_R_U ] =  R_U + i2_U_out            - r_U_to_V
        
        ## Vaccinateds:
        curr[...,  self.ix_S_V] =  S_V - s_V_out             + s_U_to_V
        curr[..., self.ix_E1_V] = E1_V +  s_V_out - e1_V_out + e1_U_to_V
        curr[..., self.ix_E2_V] = E2_V + e1_V_out - e2_V_out + e2_U_to_V
        curr[..., self.ix_I1_V] = I1_V + e2_V_out - i1_V_out + i1_U_to_V
        curr[..., self.ix_I2_V] = I2_V + i1_V_out - i2_V_out + i2_U_to_V
        curr[..., self.ix_R_V ] =  R_V + i2_V_out            + r_U_to_V

        # Keep parameters fixed.
        curr[..., self.ix_R0:] = prev[..., self.ix_R0:]

        # Record the R0(t) values for each particle.
        curr[..., self.ix_R0_val] = R0

        # Keep track of adjustment to be used when forecasting.
        curr[..., self.ix_adjustment] = adjustment

        # Keep track of what Ei/Et we are using for diagnostics later
        curr[..., self.ix_mean_Ei] = mean_vacc_Ei
        curr[..., self.ix_mean_Et] = mean_vacc_Et


    def get_R0_from_lookup(self, ctx, step_date, is_fs, R0_ix, params):
        start = ctx.component['time'].start
        forecast_with_future_R0 = False
        param_name = 'forecast_with_future_R0'
        if 'model' in params and param_name in params['model']:
            forecast_with_future_R0 = params['model'][param_name]
        if is_fs and not forecast_with_future_R0:
            # NOTE: Forecasting run, only using Reff(forecast_date).
            when = start
        else:
            when = step_date
        # Retrieve R0(t) values from the lookup table.
        R0_values = self.__R0_lookup.lookup(when)
        
        if ctx.params['model']['reorder_reffs']:
            R0_ix_mapped = self.R0_order_map[R0_ix]

            return R0_values[R0_ix_mapped]
        else:
            return R0_values[R0_ix]


    def get_external_exposure_from_lookup(self, step_date, shape):
        external = np.zeros(shape)
        if self.__external_lookup is not None:
            external_values = self.__external_lookup.lookup(step_date)
            n = len(external_values)
            if n == 1:
                external[:] = external_values[0]
            elif n == len(external):
                # NOTE: we currently assume that when there are multiple
                # external exposure trajectories, that the values will only be
                # non-zero in the forecasting period (i.e., there are no more
                # observations, so particles will not be resampled) and we can
                # simply assign the trajectories to each particle in turn.
                external[:] = external_values[:]
            else:
                raise ValueError('Invalid number of lookup values: {}'
                                 .format(n))

        return external

    def get_vaccinations_from_lookup(self, step_date, shape):
        vaccinations = np.zeros(shape)
        mean_Ei = np.zeros(shape)
        mean_Et = np.zeros(shape)
        
        
        if self.__vaccinations_lookup is not None:
            vaccinations_values = self.__vaccinations_lookup.lookup(step_date)
            n = len(vaccinations_values)
            if n == 3:
                vaccinations[:] = vaccinations_values[0]
                mean_Ei[:] = vaccinations_values[1]
                mean_Et[:] = vaccinations_values[2]
                
            else:
                raise ValueError('Invalid number of lookup values in vaccinations lookup: {}'
                                 .format(n))

        return vaccinations, mean_Ei, mean_Et

    def pre_resample(self, ctx, step_date, curr):
        """
        Change curr such that resampling is performed correctly wrt. Reff trajectories.

        :param step_date: The date and time of the current time-step.
        :param curr: The model states at the end of the observation period.
        """
        if self.__R0_lookup is None:
            print("No R0_lookup!")
            return
            
        if not ctx.params['model']['reorder_reffs']:
            return

        R0_values = self.__R0_lookup.lookup(step_date)

        self.R0_order_map = np.argsort(R0_values)

        R0_ix = np.around(curr[..., self.ix_R0_ix]).astype(int)
        curr[..., self.ix_R0_ix] = self.R0_order_map[R0_ix]
        


    def pr_inf(self, prev, curr):
        """
        Return the probability of an individual becoming infected, for any
        number of state vectors.

        :param prev: The model states at the start of the observation period.
        :param curr: The model states at the end of the observation period.
        """
        # Count the number of susceptible / exposed individuals at both ends
        # of the simulation period.
        prev_amt = np.sum(prev[..., self.ix_S_U:self.ix_I2_U] + prev[..., self.ix_S_V:self.ix_I2_V], axis=-1)
        curr_amt = np.sum(curr[..., self.ix_S_U:self.ix_I2_U] + curr[..., self.ix_S_V:self.ix_I2_V], axis=-1)
        # Avoid returning very small negative values (e.g., -1e-10).
        num_infs = np.maximum(prev_amt - curr_amt, 0)
        return num_infs / self.popn_size

    def is_seeded(self, hist):
        """Identify state vectors where infections have occurred.

        :param hist: A matrix of arbitrary dimensions, whose final dimension
            covers the model state space (i.e., has a length no smaller than
            that returned by :py:func:`state_size`).
        :type hist: numpy.ndarray

        :returns: A matrix of one fewer dimensions than ``hist`` that contains
            ``1`` for state vectors where infections have occurred and ``0``
            for state vectors where they have not.
        :rtype: numpy.ndarray
        """
        return np.ceil(1 - hist[..., self.ix_S_U + self.ix_S_V])

    def is_extinct(self, hist):
        """
        Return an array that identifies state vectors where the epidemic has
        become extinct.

        By default, this method returns ``False`` for all particles.
        Stochastic models should override this method.

        :param hist: A matrix of arbitrary dimensions, whose final dimension
            covers the model state space (i.e., has a length no smaller than
            that returned by :py:func:`state_size`).
        :type hist: numpy.ndarray

        :returns: A matrix of one fewer dimensions than ``hist`` that contains
            ``True`` for state vectors where the epidemic is extinct and
            ``False`` for state vectors where the epidemic is ongoing.
        :rtype: numpy.ndarray
        """
        # Count the number of individuals in E1, E2, I1, and I2.
        num_exposed = np.sum(hist[..., self.ix_E1_U:self.ix_R], axis=-1)
        return num_exposed == 0

    def is_valid(self, hist):
        """Ignore state vectors where no infections have occurred, as their
        properties (such as parameter distributions) are uninformative."""
        return self.is_seeded(hist)

    def describe(self):
        descr = [info_tuple for info_tuple in self.__info]
        # Check whether R0_ix can be smoothed by, e.g., post-regularistion.
        if self.__regularise_R0_ix:
            for ix in range(len(descr)):
                if descr[ix][0] == 'R0_ix':
                    descr[ix] = (descr[ix][0], True, *descr[ix][2:])
        return descr

    def stat_info(self):
        """
        Return the summary statistics that are provided by this model.

        Each statistic is represented as a ``(name, stat_fn)`` tuple, where
        ``name`` is a string and ``stat_fn`` is a function that accepts one
        argument (the particle history matrix) and returns the statistic (see,
        e.g., :py:func:`stat_generation_interval`).
        """
        return [("gen_int", self.stat_generation_interval)]

    def stat_generation_interval(self, hist):
        """
        Calculate the mean generation interval for each particle.

        :param hist: The particle history matrix, or a subset thereof.
        """
        return 1 / hist[..., self.ix_sigma] + 0.75 / hist[..., self.ix_gamma]
