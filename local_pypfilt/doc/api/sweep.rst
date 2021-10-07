pypfilt.sweep
=============

.. py:module:: pypfilt.sweep

The :mod:`pypfilt.sweep` module iterates over forecast scenarios defined in
any number of scenario files.

Forecast scenarios
------------------

.. autofunction:: forecasts

.. autoclass:: Forecasts

Forecasting in parallel
-----------------------

The :mod:`pypfilt.sweep` module also provides support for running forecasts
across multiple Python processes:

.. autofunction:: forecasts_mp

.. autofunction:: get_forecasts_mp
