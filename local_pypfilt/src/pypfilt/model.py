"""Base class for simulation models."""

import abc
import numpy as np


class Model(abc.ABC):
    """
    The base class for simulation models, which defines the minimal set of
    methods that are required.
    """

    @abc.abstractmethod
    def init(self, ctx, vec):
        """
        Initialise a matrix of state vectors.

        :param ctx: The simulation context.
        :param vec: An uninitialised :math:`P \\times S` matrix of state
            vectors, for :math:`P` particles and state vectors of length
            :math:`S` (as defined by :py:func:`~state_size`).
            To set, e.g., the first element of each state vector to :math:`1`,
            you can use an ellipsis slice: :code:`vec[..., 0] = 1`.
        """
        pass

    @abc.abstractmethod
    def state_size(self):
        """
        Return the size of the state vector.
        """
        pass

    @abc.abstractmethod
    def update(self, params, step_date, dt, is_fs, prev, curr):
        """
        Perform a single time-step.

        :param params: Simulation parameters.
        :param step_date: The date and time of the current time-step.
        :param dt: The time-step size (days).
        :param is_fs: Indicates whether this is a forecasting simulation.
        :param prev: The state before the time-step.
        :param curr: The state after the time-step (destructively updated).
        """
        pass

    @abc.abstractmethod
    def pre_resample(self, ctx, step_date, curr):
        """
        Performs pre-resampling steps.

        :param step_date: The date and time of the current time-step.
        :param curr: The state after the time-step (destructively updated).
        """
        pass

    @abc.abstractmethod
    def describe(self):
        """
        Describe each component of the state vector with a tuple of the form
        ``(name, smooth, min, max)``, where ``name`` is a descriptive name for
        the variable/parameter, ``smooth`` is a boolean that indicates whether
        the parameter admits continuous sampling (e.g., post-regularisation),
        and ``min`` and ``max`` define the (inclusive) range of valid values.
        These tuples **must** be in the same order as the state vector itself.
        """

    def resume_from_cache(self, ctx):
        """
        Notify the model that a simulation will begin from a saved state.

        The model does not need to initialise the state vectors, since these
        will have been loaded from a cache file, but it may need to update any
        internal variables (i.e., those not stored in the state vectors).

        .. note:: Models should only implement this method if they need to
           prepare for the simulation.
        """
        pass

    def stat_info(self):
        """
        Describe each statistic that can be calculated by this model as a
        ``(name, stat_fn)`` tuple, where ``name`` is a string that identifies
        the statistic and ``stat_fn`` is a function that calculates the value
        of the statistic.

        .. note:: Models should only implement this method if they define one
           or more statistics.
        """
        return []

    def is_valid(self, hist):
        """
        Identify particles whose state and parameters can be inspected. By
        default, this function returns ``True`` for all particles. Override
        this function to ensure that inchoate particles are correctly
        ignored.

        .. note:: Models should only implement this method if there are
           conditions where some particles should be ignored.
        """
        return np.ones((hist.shape[-2],), dtype=bool)
