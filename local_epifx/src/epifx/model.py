"""Models of disease transmission in human populations."""

import abc
import numpy as np
import pypfilt


class Model(pypfilt.model.Model):
    """
    The base class for simulation models that are used with observation models
    and summary tables provided by ``epifx``.
    """

    @abc.abstractmethod
    def sample_columns(self):
        """
        Identify the parameters in the model state vector that can be saved
        to, and loaded from, external data files.

        This should be returned as a dictionary that maps parameter names
        (which must not contain spaces) to state vector indices.
        """
        pass

    @abc.abstractmethod
    def population_size(self):
        """
        Return the model population size.

        This is used by ``epifx.obs.PopnCounts`` to calculate expected
        observation values.
        """
        pass

    @abc.abstractmethod
    def pr_inf(self, prev, curr):
        """
        Return the likelihood of an individual becoming infected, for
        number of state vectors.

        :param prev: The model states at the start of the observation period.
        :param curr: The model states at the end of the observation period.

        This is used by ``epifx.obs.PopnCounts`` and
        ``epifx.obs.SampleCounts`` to calculate expected observation values.
        """
        pass

    @abc.abstractmethod
    def is_seeded(self, hist):
        """
        Return an array that identifies state vectors where infections have
        occurred.

        :param hist: A matrix of arbitrary dimensions, whose final dimension
            covers the model state space (i.e., has a length no smaller than
            that returned by :py:func:`state_size`).
        :type hist: numpy.ndarray

        :returns: A matrix of one fewer dimensions than ``hist`` that contains
            ``1`` for state vectors where infections have occurred and ``0``
            for state vectors where they have not.
        :rtype: numpy.ndarray

        This is used by ``epifx.summary.PrOutbreak`` to calculate outbreak
        probabilities.
        """
        pass

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
        shape = hist[..., 0].shape
        return np.full(shape, False)
