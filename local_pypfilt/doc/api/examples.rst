pypfilt.examples
================

.. py:module:: pypfilt.examples

The :mod:`pypfilt.examples` module provides an example
:ref:`predator-prey system <gs-equations>`.

Predator-prey system
--------------------

.. py:module:: pypfilt.examples.predation

Models
^^^^^^

.. autoclass:: LotkaVolterra
   :members:

.. autoclass:: ObsModel
   :members:

Example files
^^^^^^^^^^^^^

.. autofunction:: write_example_files

.. autofunction:: example_toml_data

.. autofunction:: example_obs_x_data

.. autofunction:: example_obs_y_data

Generating forecasts
^^^^^^^^^^^^^^^^^^^^

.. autofunction:: forecast

.. autofunction:: plot

.. autofunction:: plot_params

.. autofunction:: plot_forecasts


Other functions
^^^^^^^^^^^^^^^^^^^^

.. autofunction:: default_priors

.. autofunction:: make_params

.. autofunction:: save_scalar_observations
