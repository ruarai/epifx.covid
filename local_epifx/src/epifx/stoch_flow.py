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
        ("R_bias", False, 1/3, 2.0),
        ("adjustment", False, 0, 1)]

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

    ix_R_bias = 18
    ix_adjustment = 19

    R0_order_map = np.arange(0, 1000, 1)

    sigma_transitions = None
    gamma_transitions = None
    vacc_transitions = None

    comp_mask_all = None
    comp_mask_U = None

    
    n_compartments = 12


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


        sigma_transitions_ix = [(self.ix_E1_U, self.ix_E2_U), (self.ix_E2_U, self.ix_I1_U),
                                (self.ix_E1_V, self.ix_E2_V), (self.ix_E2_V, self.ix_I1_V)]
        self.sigma_transitions = np.moveaxis(np.array(sigma_transitions_ix), -1, 0)

        gamma_transitions_ix = [(self.ix_I1_U, self.ix_I2_U), (self.ix_I2_U, self.ix_R_U),
                                (self.ix_I1_V, self.ix_I2_V), (self.ix_I2_V, self.ix_R_V)]
        self.gamma_transitions = np.moveaxis(np.array(gamma_transitions_ix), -1, 0)

        vacc_transitions_ix = [(self.ix_S_U, self.ix_S_V), (self.ix_E1_U, self.ix_E1_V), (self.ix_E2_U, self.ix_E2_V),
                               (self.ix_I1_U, self.ix_I1_V), (self.ix_I2_U, self.ix_I2_V), (self.ix_R_U, self.ix_R_V)]
        self.vacc_transitions = np.moveaxis(np.array(vacc_transitions_ix), -1, 0)

                
        self.comp_mask_U = np.array([1,1,1,1,1,1,0,0,0,0,0,0])
        self.comp_mask_all = np.ones(self.n_compartments)




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
            'R_bias': self.ix_R_bias,
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

        R0_ix = np.around(prev[..., self.ix_R0_ix]).astype(int)
        R_bias = prev[..., self.ix_R_bias].copy()
        adjustment = prev[..., self.ix_adjustment].copy()

        current_state = prev[:, 0:(self.ix_R_V + 1)]

        
        epoch = ctx.component['time'].to_scalar(ctx.params['epoch'])
        curr_t = ctx.component['time'].to_scalar(step_date)
        zero_mask = prev[..., self.ix_t0] > (curr_t - epoch)

        if self.__R0_lookup is not None:
            R0 = self.get_R0_from_lookup(ctx, step_date, is_fs, R0_ix, params)

        external_exp = self.get_external_exposure_from_lookup(step_date, R0.shape)
        vacc_rate, mean_vacc_Ei, mean_vacc_Et = self.get_vaccinations_from_lookup(step_date, R0.shape)

        n_pop = self.popn_size


        # Only calculate adj. factor for backcasts (expecting Reff to be held constant across forecasting period)
        if not is_fs:

            I_U, I_V = prev[..., self.ix_I1_U] + prev[..., self.ix_I2_U], prev[..., self.ix_I1_V] + prev[..., self.ix_I2_V]
            S_U, S_V = prev[..., self.ix_S_U], prev[..., self.ix_S_V]
                
            # Calculating our adjustment factor
            # Not factoring out the 1/n out for numerical stability
            adjustment = (I_U / n_pop + I_V / n_pop) / ((I_U / n_pop + (1 - mean_vacc_Et) * I_V / n_pop) * (S_U / n_pop + (1 - mean_vacc_Ei) * S_V / n_pop))
            adjustment = np.nan_to_num(adjustment, nan = 0, posinf = 0)

        R0[zero_mask] = 0
        sigma[zero_mask] = 0
        gamma[zero_mask] = 0

        beta = R0 * adjustment * np.exp2(R_bias) * gamma


        n_U = np.dot(current_state, self.comp_mask_U)
        vacc_increase = np.divide(vacc_rate, n_U, out=np.zeros_like(n_U), where = n_U != 0)


        compartment_out, compartment_in = self.get_compartment_change(current_state,
                                                                      sigma, gamma, beta, external_exp, vacc_increase, 
                                                                      mean_vacc_Et, mean_vacc_Ei, n_pop,
                                                                      rnd, dt)


        curr[..., 0:(self.ix_R_V + 1)] = current_state + compartment_in - compartment_out # Flow between compartments
        curr[..., self.ix_R0:] = prev[..., self.ix_R0:] # Keep parameters fixed.
        curr[..., self.ix_R0_val] = R0 # Record the R0(t) values for each particle.
        curr[..., self.ix_adjustment] = adjustment # Keep track of adjustment to be used when forecasting.

    def get_compartment_change(self, current_state, sigma, gamma, beta, external_exp, vacc_increase, mean_vacc_Et, mean_vacc_Ei, n_pop, rnd, dt):

        lambda_inf = current_state[..., self.ix_I1_U] + current_state[..., self.ix_I2_U] + \
                            (1 - mean_vacc_Et) * (current_state[..., self.ix_I1_V] + current_state[..., self.ix_I2_V])

        n_particles = sigma.shape[0]

        # flow_rate defines the transition rates across each n_particle particle, from one compartment to another
        flow_rate = np.zeros((n_particles , self.n_compartments, self.n_compartments))

        flow_rate[:, self.ix_S_U,  self.ix_E1_U] = (beta * lambda_inf + external_exp) / n_pop
        flow_rate[:, self.ix_S_V,  self.ix_E1_V] = beta * lambda_inf * (1 - mean_vacc_Ei) / n_pop

        flow_rate[np.index_exp[:] + tuple(self.sigma_transitions)] = 2 * sigma[:,None]
        flow_rate[np.index_exp[:] + tuple(self.gamma_transitions)] = 2 * gamma[:,None]
        flow_rate[np.index_exp[:] + tuple(self.vacc_transitions)]  = vacc_increase[:,None]



        flow_rate = -np.expm1(-flow_rate * dt)
        flow_rate = rnd.binomial(current_state[..., None].astype(int), flow_rate)

        # Scale down our outflow if > our compartment counts
        out_row_sums = np.sum(flow_rate, axis = 2)
        excess_out = (out_row_sums > current_state)

        # Calculate outflows as proportion of max. possible
        as_proportions = flow_rate[excess_out, :].astype(float) / out_row_sums[excess_out, None]

        # Scale down to fill max possible
        flow_rate[excess_out, :] = np.around(as_proportions * current_state[excess_out, None]).astype(int)


        compartment_out = np.sum(flow_rate, axis = 2) # Replace with np.dot sometime
        compartment_in = np.sum(flow_rate, axis = 1)

        return compartment_out, compartment_in


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
