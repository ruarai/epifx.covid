API documentation
=================

The :mod:`pypfilt` module provides top-level functions for running forecasts
and simulating observations from simulation models:

:func:`~pypfilt.forecasts_iter`
   Iterate over forecast settings for each scenario.

:func:`~pypfilt.forecast`
   Generate forecasts at specific times for a single scenario.

:func:`~pypfilt.fit`
   Fit the simulation model to all of the available observations.

:func:`~pypfilt.simulate_from_model`
   Simulate observations from the simulation model, according to each
   observation model.

It also contains a number of sub-modules.
Some are intended for public use (see the :ref:`key modules <tbl_key_mods>`
table), while others are likely of no use outside of pypfilt_ (see the
:ref:`secondary modules <tbl_other_mods>` table).

.. table:: Key pypfilt_ modules.
   :widths: auto
   :name: tbl_key_mods

   =======================  ==================================================
   Module                   Description
   =======================  ==================================================
   :mod:`pypfilt`           Provides model-fitting and forecasting functions
   :mod:`pypfilt.config`    Reads forecast scenarios from TOML_ files
   :mod:`pypfilt.sweep`     Iterates over forecast scenarios
   :mod:`pypfilt.model`     Defines the simulation model base class
                            :class:`~pypfilt.model.Model`
   :mod:`pypfilt.obs`       Defines the observation model base class
                            :class:`~pypfilt.obs.Obs`
   :mod:`pypfilt.params`    Provides default parameter values for the particle
                            filter
   :mod:`pypfilt.time`      Provides scalar and date-time simulation time
                            scales
   :mod:`pypfilt.summary`   Provides common summary statistics and records
                            outputs
   :mod:`pypfilt.plot`      Provides functions for plotting summary statistics
   :mod:`pypfilt.io`        Reads data tables from text files
   :mod:`pypfilt.examples`  Provides example models
   =======================  ==================================================

.. table:: Secondary pypfilt_ modules, which you will rarely (if ever) use
           directly.
   :widths: auto
   :name: tbl_other_mods

   =======================  ==================================================
   Module                   Description
   =======================  ==================================================
   :mod:`pypfilt.cache`     Implements the particle filter state cache
   :mod:`pypfilt.check`     Provides invariant checks for the state history
                            matrix
   :mod:`pypfilt.context`   Manages simulation components, parameters, event
                            handlers, etc.
   :mod:`pypfilt.pfilter`   The particle filter core: time-steps and adjusting
                            particle weights
   :mod:`pypfilt.resample`  Implements particle resampling and
                            post-regularisation
   :mod:`pypfilt.state`     Creates the state history matrix
   :mod:`pypfilt.stats`     Calculates weighted quantiles, credible intervals,
                            etc
   =======================  ==================================================

.. toctree::
   :maxdepth: 2
   :caption: Key modules

   pypfilt
   config
   sweep
   model
   obs
   params
   time
   summary
   plot
   io
   examples

.. toctree::
   :maxdepth: 2
   :caption: Secondary modules

   cache
   check
   context
   pfilter
   resample
   state
   stats
