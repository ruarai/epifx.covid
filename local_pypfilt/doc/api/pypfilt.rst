pypfilt
=======

.. py:module:: pypfilt

The :mod:`pypfilt` module provides top-level functions for running forecasts
and simulating observations from simulation models.

.. autofunction:: pypfilt.forecasts_iter

.. autofunction:: pypfilt.forecast

.. autofunction:: pypfilt.fit

.. autofunction:: pypfilt.simulate_from_model

The :mod:`pypfilt` module also re-exports a number of items from sub-modules:

* Model base classes: :class:`~pypfilt.model.Model` and
  :class:`~pypfilt.obs.Obs`.

* Summary statistic base classes: :class:`~pypfilt.summary.Monitor` and
  :class:`~pypfilt.summary.Table`.

* Simulation time scales: :class:`~pypfilt.time.Datetime` and
  :class:`~pypfilt.time.Scalar`.

* Configuration objects: :class:`~pypfilt.config.Config` and
  :class:`~pypfilt.context.Context`.

* Default particle filter parameters: :func:`~pypfilt.params.default_params`.
