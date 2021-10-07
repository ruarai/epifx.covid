"""Deterministic models of infectious diseases."""

import logging
from pathlib import Path
import numpy as np
import pypfilt
from .model import Model


class SEIR(Model):
    r"""
    An SEIR compartment model for a single circulating influenza strain, under
    the assumption that recovered individuals are completely protected against
    reinfection.

    .. math::

        \frac{dS}{dt} &= - \beta S^\eta I \\[0.5em]
        \frac{dE}{dt} &= \beta S^\eta I - \sigma E \\[0.5em]
        \frac{dI}{dt} &= \sigma E - \gamma I \\[0.5em]
        \frac{dR}{dt} &= \gamma I \\[0.5em]
        \beta &= R_0 \cdot \gamma

    ==============  ================================================
    Parameter       Meaning
    ==============  ================================================
    :math:`R_0`     Basic reproduction number
    :math:`\sigma`  Inverse of the incubation period (day :sup:`-1`)
    :math:`\gamma`  Inverse of the infectious period (day :sup:`-1`)
    :math:`\eta`    Inhomogeneous social mixing coefficient
    :math:`\alpha`  Temporal forcing coefficient
    ==============  ================================================

    The force of infection can be subject to temporal forcing :math:`F(t)`, as
    mediated by :math:`\alpha`:

    .. math::

        \beta(t) = \beta \cdot \left[1 + \alpha \cdot F(t)\right]

    Note that this requires the forcing time-series to be stored in the lookup
    table ``'R0_forcing'``.
    """

    __info = [("S", False, 0, 1), ("E", False, 0, 1), ("I", False, 0, 1),
              ("R", False, 0, 1),
              ("R0", True, 1, 2), ("sigma", True, 1/3, 2),
              ("gamma", True, 1/3, 1), ("eta", True, 1, 2),
              ("alpha", True, -0.2, 0.2),
              ("t0", False, 0, 50),
              ]

    ix_S = 0
    ix_E = 1
    ix_I = 2
    ix_R = 3
    ix_R0 = 4
    ix_sigma = 5
    ix_gamma = 6
    ix_eta = 7
    ix_alpha = 8
    ix_t0 = 9

    def __init__(self):
        """Initialise the model instance."""
        self.__Forcing_lookup = None

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
        self.__Forcing_lookup = None

        prior = ctx.params['model']['prior']
        rnd = ctx.component['random']['model']
        rnd_size = vec[..., 0].shape

        # Initialise the model state (fully susceptible population).
        initial_exposures = 1.0 / self.popn_size
        vec[..., :] = 0
        vec[..., self.ix_S] = 1 - initial_exposures
        vec[..., self.ix_E] = initial_exposures
        vec[..., self.ix_R0] = prior['R0'](rnd, size=rnd_size)
        vec[..., self.ix_sigma] = prior['sigma'](rnd, size=rnd_size)
        vec[..., self.ix_gamma] = prior['gamma'](rnd, size=rnd_size)
        vec[..., self.ix_eta] = prior['eta'](rnd, size=rnd_size)
        vec[..., self.ix_alpha] = 0
        vec[..., self.ix_t0] = prior['t0'](rnd, size=rnd_size)

        self.load_samples_file(ctx, vec)
        self.load_lookup_tables(ctx, vec, init_values=True)

    def sample_columns(self):
        """Identify parameters that can be saved and loaded."""
        ix_tbl = {
            'R0': self.ix_R0,
            'sigma': self.ix_sigma,
            'gamma': self.ix_gamma,
            'eta': self.ix_eta,
            'alpha': self.ix_alpha,
            't0': self.ix_t0,
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

    def load_lookup_tables(self, ctx, vec, init_values=False):
        logger = logging.getLogger(__name__)
        rnd = ctx.component['random']['model']
        rnd_size = vec[..., 0].shape
        prior = ctx.params['model']['prior']
        tables = ctx.component.get('lookup', {})
        if 'R0_forcing' in tables and self.__Forcing_lookup is None:
            self.__Forcing_lookup = tables['R0_forcing']
            logger.info('Using lookup table for R0 forcing with {} values'
                        .format(self.__Forcing_lookup.value_count()))
            if init_values:
                vec[..., self.ix_alpha] = prior['alpha'](rnd, size=rnd_size)

    def update(self, ctx, step_date, dt, is_fs, prev, curr):
        """Perform a single time-step.

        :param ctx: The simulation context.
        :param step_date: The date and time of the current time-step.
        :param dt: The time-step size (days).
        :param is_fs: Indicates whether this is a forecasting simulation.
        :param prev: The state before the time-step.
        :param curr: The state after the time-step (destructively updated).
        """
        # Update parameters and lookup tables that are defined in self.init()
        # and which will not exist if we are resuming from a cached state.
        self.popn_size = ctx.params['model']['population_size']
        self.load_lookup_tables(ctx, prev, init_values=False)

        # Extract each parameter.
        R0 = prev[..., self.ix_R0].copy()
        sigma = prev[..., self.ix_sigma].copy()
        gamma = prev[..., self.ix_gamma].copy()
        eta = prev[..., self.ix_eta].copy()
        alpha = prev[..., self.ix_alpha].copy()
        t0 = prev[..., self.ix_t0].copy()

        beta = R0 * gamma
        if self.__Forcing_lookup is not None:
            # Modulate the force of infection with temporal forcing.
            force = alpha * self.__Forcing_lookup.lookup(step_date)[0]
            # Ensure the force of infection is non-negative (can be zero).
            beta *= np.maximum(1.0 + force, 0)

        epoch = ctx.component['time'].to_scalar(ctx.params['epoch'])
        curr_t = ctx.component['time'].to_scalar(step_date)
        zero_mask = t0 > (curr_t - epoch)
        R0[zero_mask] = 0
        sigma[zero_mask] = 0
        gamma[zero_mask] = 0
        eta[zero_mask] = 0
        alpha[zero_mask] = 0
        t0[zero_mask] = 0
        beta[zero_mask] = 0

        # Extract each compartment.
        S = prev[..., self.ix_S]
        E = prev[..., self.ix_E]
        I = prev[..., self.ix_I]
        R = prev[..., self.ix_R]

        # Calculate flows between compartments.
        s_to_e = dt * beta * I * S ** eta
        e_to_i = dt * sigma * E
        i_to_r = dt * gamma * I

        # Update the compartment values.
        curr[..., self.ix_S] = S - s_to_e
        curr[..., self.ix_E] = E + s_to_e - e_to_i
        curr[..., self.ix_I] = I + e_to_i - i_to_r

        # Enforce invariants on the S, E, and I compartments.
        curr[..., :self.ix_R] = np.clip(curr[..., :self.ix_R], 0, 1)
        sum_SEI = np.sum(curr[..., :self.ix_R], axis=-1)
        mask_invalid = sum_SEI > 1
        if np.any(mask_invalid):
            denom = sum_SEI[mask_invalid, None]
            curr[mask_invalid, :self.ix_R] = (curr[mask_invalid, :self.ix_R]
                                              / denom)

        # Calculate the size of the R compartment and clip appropriately.
        curr[..., self.ix_R] = np.clip(
            1.0 - np.sum(curr[..., :self.ix_R], axis=-1),
            0.0, 1.0)

        # Keep parameters fixed.
        curr[..., self.ix_R0:] = prev[..., self.ix_R0:]

    def pr_inf(self, prev, curr):
        """Calculate the likelihood of an individual becoming infected, for
        any number of state vectors.

        :param prev: The model states at the start of the observation period.
        :param curr: The model states at the end of the observation period.
        """
        # Count the number of susceptible / exposed individuals at both ends
        # of the simulation period.
        prev_amt = np.sum(prev[..., 0:2], axis=-1)
        curr_amt = np.sum(curr[..., 0:2], axis=-1)
        # Avoid returning very small negative values (e.g., -1e-10).
        return np.maximum(prev_amt - curr_amt, 0)

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
        return np.ceil(1 - hist[..., 0])

    def is_valid(self, hist):
        """Ignore state vectors where no infections have occurred, as their
        properties (such as parameter distributions) are uninformative."""
        return self.is_seeded(hist)

    def describe(self):
        return self.__info

    def stat_info(self):
        """Return the details of each statistic that can be calculated by this
        model. Each such statistic is represented as a ``(name, stat_fn)``
        pair, where ``name`` is a string that identifies the statistic and
        ``stat_fn`` is a function that calculates the statistic (see, e.g.,
        :py:func:`stat_Reff`).
        """
        return [("Reff", self.stat_Reff)]

    def stat_Reff(self, hist):
        """Calculate the effective reproduction number :math:`R_{eff}` for
        every particle.

        :param hist: The particle history matrix, or a subset thereof.
        """
        return hist[..., self.ix_S] * hist[..., self.ix_R0]


class SEEIIR(Model):
    r"""An SEEIIR compartment model for a single circulating influenza strain,
    under the assumption that recovered individuals are completely protected
    against reinfection.

    .. math::

        \frac{dS}{dt} &= - \beta S^\eta (I_1 + I_2) \\[0.5em]
        \frac{dE_1}{dt} &= \beta S^\eta (I_1 + I_2) - 2 \sigma E_1 \\[0.5em]
        \frac{dE_2}{dt} &= 2 \sigma E_1 - 2 \sigma E_2 \\[0.5em]
        \frac{dI_1}{dt} &= 2 \sigma E_2 - 2 \gamma I_1 \\[0.5em]
        \frac{dI_2}{dt} &= 2 \gamma I_1 - 2 \gamma I_2 \\[0.5em]
        \frac{dR}{dt} &= 2 \gamma I_2 \\[0.5em]
        \beta &= R_0 \cdot \gamma

    ==============  ================================================
    Parameter       Meaning
    ==============  ================================================
    :math:`R_0`     Basic reproduction number
    :math:`\sigma`  Inverse of the incubation period (day :sup:`-1`)
    :math:`\gamma`  Inverse of the infectious period (day :sup:`-1`)
    :math:`\eta`    Inhomogeneous social mixing coefficient
    :math:`\alpha`  Temporal forcing coefficient
    ==============  ================================================

    The force of infection can be subject to temporal forcing :math:`F(t)`, as
    mediated by :math:`\alpha`:

    .. math::

        \beta(t) = \beta \cdot \left[1 + \alpha \cdot F(t)\right]

    Note that this requires the forcing time-series to be stored in the lookup
    table ``'R0_forcing'``.
    """

    __info = [
        ("S", False, 0, 1), ("E1", False, 0, 1), ("E2", False, 0, 1),
        ("I1", False, 0, 1), ("I2", False, 0, 1), ("R", False, 0, 1),
        ("R0", True, 1, 2), ("sigma", True, 1/3, 2),
        ("gamma", True, 1/3, 1), ("eta", True, 1, 2),
        ("alpha", True, -0.2, 0.2),
        ("t0", False, 0, 100)]

    ix_S = 0
    ix_E1 = 1
    ix_E2 = 2
    ix_I1 = 3
    ix_I2 = 4
    ix_R = 5
    ix_R0 = 6
    ix_sigma = 7
    ix_gamma = 8
    ix_eta = 9
    ix_alpha = 10
    ix_t0 = 11

    def __init__(self):
        """Initialise the model instance."""
        self.__R0_lookup = None
        self.__Overseas_lookup = None
        self.__Forcing_lookup = None

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
        self.__Overseas_lookup = None
        self.__Forcing_lookup = None

        prior = ctx.params['model']['prior']
        rnd = ctx.component['random']['model']
        rnd_size = vec[..., 0].shape

        # Initialise the model state (fully susceptible population).
        initial_exposures = 1.0 / self.popn_size
        vec[..., :] = 0
        vec[..., self.ix_S] = 1 - initial_exposures
        vec[..., self.ix_E1] = initial_exposures
        vec[..., self.ix_R0] = prior['R0'](rnd, size=rnd_size)
        vec[..., self.ix_sigma] = prior['sigma'](rnd, size=rnd_size)
        vec[..., self.ix_gamma] = prior['gamma'](rnd, size=rnd_size)
        vec[..., self.ix_eta] = prior['eta'](rnd, size=rnd_size)
        vec[..., self.ix_alpha] = 0
        vec[..., self.ix_t0] = prior['t0'](rnd, size=rnd_size)

        self.load_samples_file(ctx, vec)
        self.load_lookup_tables(ctx, vec, init_values=True)

    def sample_columns(self):
        """Identify parameters that can be saved and loaded."""
        ix_tbl = {
            'R0': self.ix_R0,
            'sigma': self.ix_sigma,
            'gamma': self.ix_gamma,
            'eta': self.ix_eta,
            'alpha': self.ix_alpha,
            't0': self.ix_t0,
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

    def load_lookup_tables(self, ctx, vec, init_values=False):
        logger = logging.getLogger(__name__)
        rnd = ctx.component['random']['model']
        rnd_size = vec[..., 0].shape
        prior = ctx.params['model']['prior']
        tables = ctx.component.get('lookup', {})
        if 'R0' in tables and self.__R0_lookup is None:
            # TODO: R0_ix and R0_val
            self.__R0_lookup = tables['R0']
            num_values = self.__R0_lookup.value_count()
            logger.info('Using lookup table for R0 with {} values'.format(
                num_values))
            if init_values and num_values > 1:
                vec[..., self.ix_R0_ix] = rnd.integers(num_values,
                                                       size=rnd_size)
            elif init_values:
                vec[..., self.ix_R0_ix] = 0
        if 'Overseas Cases' in tables and self.__Overseas_lookup is None:
            self.__Overseas_lookup = tables['Overseas Cases']
            logger.info('Using lookup table for overseas cases with {} values'
                        .format(self.__Overseas_lookup.value_count()))
        if 'R0_forcing' in tables and self.__Forcing_lookup is None:
            self.__Forcing_lookup = tables['R0_forcing']
            logger.info('Using lookup table for R0 forcing with {} values'
                        .format(self.__Forcing_lookup.value_count()))
            if init_values:
                vec[..., self.ix_alpha] = prior['alpha'](rnd, size=rnd_size)

    def update(self, ctx, step_date, dt, is_fs, prev, curr):
        """Perform a single time-step.

        :param ctx: The simulation context.
        :param step_date: The date and time of the current time-step.
        :param dt: The time-step size (days).
        :param is_fs: Indicates whether this is a forecasting simulation.
        :param prev: The state before the time-step.
        :param curr: The state after the time-step (destructively updated).
        """
        # Update parameters and lookup tables that are defined in self.init()
        # and which will not exist if we are resuming from a cached state.
        self.popn_size = ctx.params['model']['population_size']
        self.load_lookup_tables(ctx, prev, init_values=False)

        # Extract each parameter.
        R0 = prev[..., self.ix_R0].copy()
        sigma = prev[..., self.ix_sigma].copy()
        gamma = prev[..., self.ix_gamma].copy()
        eta = prev[..., self.ix_eta].copy()
        alpha = prev[..., self.ix_alpha].copy()
        t0 = prev[..., self.ix_t0].copy()

        # if self.__R0_lookup is not None:
        #     start = ctx.component['time'].start
        #     forecast_with_future_R0 = False
        #     param_name = 'forecast_with_future_R0'
        #     if 'model' in params and param_name in params['model']:
        #         forecast_with_future_R0 = params['model'][param_name]
        #     if is_fs and not forecast_with_future_R0:
        #         # NOTE: Forecasting run, only using Reff(forecast_date).
        #         when = start
        #     else:
        #         when = step_date
        #     # Retrieve R0(t) values from the lookup table.
        #     R0_values = self.__R0_lookup.lookup(when)
        #     R0 = R0_values[R0_ix]

        beta = R0 * gamma
        if self.__Forcing_lookup is not None:
            # Modulate the force of infection with temporal forcing.
            force = alpha * self.__Forcing_lookup.lookup(step_date)[0]
            # Ensure the force of infection is non-negative (can be zero).
            beta *= np.maximum(1.0 + force, 0)

        import_rate = 0
        if self.__Overseas_lookup is not None:
            imports = self.__Overseas_lookup.lookup(step_date)[0]
            import_rate = imports / self.popn_size

        epoch = ctx.component['time'].to_scalar(ctx.params['epoch'])
        curr_t = ctx.component['time'].to_scalar(step_date)
        zero_mask = t0 > (curr_t - epoch)
        R0[zero_mask] = 0
        sigma[zero_mask] = 0
        gamma[zero_mask] = 0
        eta[zero_mask] = 0
        alpha[zero_mask] = 0
        t0[zero_mask] = 0
        beta[zero_mask] = 0

        # Extract each compartment.
        S = prev[..., self.ix_S]
        E1 = prev[..., self.ix_E1]
        E2 = prev[..., self.ix_E2]
        I1 = prev[..., self.ix_I1]
        I2 = prev[..., self.ix_I2]

        # Calculate flows between compartments.
        s_to_e1 = dt * (beta * (I1 + I2) * S ** eta + import_rate)
        e1_to_e2 = dt * 2 * sigma * E1
        e2_to_i1 = dt * 2 * sigma * E2
        i1_to_i2 = dt * 2 * gamma * I1
        i2_to_r = dt * 2 * gamma * I2

        # Update the compartment values.
        curr[..., self.ix_S] = S - s_to_e1
        curr[..., self.ix_E1] = E1 + s_to_e1 - e1_to_e2
        curr[..., self.ix_E2] = E2 + e1_to_e2 - e2_to_i1
        curr[..., self.ix_I1] = I1 + e2_to_i1 - i1_to_i2
        curr[..., self.ix_I2] = I2 + i1_to_i2 - i2_to_r

        # Enforce invariants on the S, E, and I compartments.
        curr[..., :self.ix_R] = np.clip(curr[..., :self.ix_R], 0, 1)
        sum_SEI = np.sum(curr[..., :self.ix_R], axis=-1)
        mask_invalid = sum_SEI > 1
        if np.any(mask_invalid):
            denom = sum_SEI[mask_invalid, None]
            curr[mask_invalid, :self.ix_R] = (curr[mask_invalid, :self.ix_R]
                                              / denom)

        # Calculate the size of the R compartment and clip appropriately.
        curr[..., self.ix_R] = np.clip(
            1.0 - np.sum(curr[..., :self.ix_R], axis=-1),
            0.0, 1.0)

        # Keep parameters fixed.
        curr[..., self.ix_R0:] = prev[..., self.ix_R0:]

    def pr_inf(self, prev, curr):
        """Calculate the likelihood of an individual becoming infected, for
        any number of state vectors.

        :param prev: The model states at the start of the observation period.
        :param curr: The model states at the end of the observation period.
        """
        # Count the number of susceptible / exposed individuals at both ends
        # of the simulation period.
        prev_amt = np.sum(prev[..., :self.ix_I2], axis=-1)
        curr_amt = np.sum(curr[..., :self.ix_I2], axis=-1)
        # Avoid returning very small negative values (e.g., -1e-10).
        return np.maximum(prev_amt - curr_amt, 0)

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
        return np.ceil(1 - hist[..., self.ix_S])

    def is_valid(self, hist):
        """Ignore state vectors where no infections have occurred, as their
        properties (such as parameter distributions) are uninformative."""
        return self.is_seeded(hist)

    def describe(self):
        return self.__info

    def stat_info(self):
        """Return the details of each statistic that can be calculated by this
        model. Each such statistic is represented as a ``(name, stat_fn)``
        pair, where ``name`` is a string that identifies the statistic and
        ``stat_fn`` is a function that calculates the statistic (see, e.g.,
        :py:func:`stat_Reff`).
        """
        return [("Reff", self.stat_Reff)]

    def stat_Reff(self, hist):
        """Calculate the effective reproduction number :math:`R_{eff}` for
        every particle.

        :param hist: The particle history matrix, or a subset thereof.
        """
        return hist[..., self.ix_S] * hist[..., self.ix_R0]
