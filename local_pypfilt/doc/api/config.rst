pypfilt.config
==============

.. py:module:: pypfilt.config

The :mod:`pypfilt.config` module allows users to define forecast scenarios in
plain-text TOML_ files.

Reading forecast scenarios
--------------------------

You can use either of the following functions to read forecast scenarios.

.. autofunction:: pypfilt.config.from_file

.. autofunction:: pypfilt.config.from_string

.. note:: You should never need to directly interact with :class:`Config`
          values returned by these two functions.
          The :mod:`pypfilt` and :mod:`pypfilt.sweep` modules provide
          convenience functions for using them.

Forecast scenario data types
----------------------------

.. autoclass:: Config

.. autoclass:: Scenario

.. autoclass:: ObsModel

.. autodata:: ParamsFn

.. autodata:: SummaryFn
