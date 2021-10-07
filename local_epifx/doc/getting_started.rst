Getting Started
===============

This page assumes that you have already :ref:`installed <install>` the
``epifx`` package, and shows how to generate forecasts using an SEIR model.

Infection models
----------------

.. autoclass:: epifx.det.SEIR
   :no-members:
   :members: init, update

.. autoclass:: epifx.det.SEEIIR
   :no-members:
   :members: init, update

.. autoclass:: epifx.stoch.SEEIIR
   :no-members:
   :members: init, update, stat_info, stat_generation_interval

.. _select:

Arbitrary model priors
----------------------

In addition to sampling each model parameter independently, the
``epifx.select`` module provides support for sampling particles according to
arbitrary target distributions, using an accept-reject sampler.

.. autofunction:: epifx.select.select

.. code-block:: python

   # Save the accepted particles to disk.
   vec = epifx.select.select(params, start, end, proposal, target, seed)

   # Load the accepted particles for use in a simulation.
   model = epifx.SEIR()
   model.set_params(vec)
   params = epifx.default_params(px_count, model, time, popn_size)


Any proposal distribution can be used with this sampler, including the default
model prior:

.. autoclass:: epifx.select.Proposal

.. autoclass:: epifx.select.DefaultProposal

Any target distribution for which a probability density can be defined can be
used with this sampler:

.. autoclass:: epifx.select.Target

Two target distributions are provided by this module.

The :class:`~TargetAny` distribution accepts all particles with equal
likelihood, for the case where the proposal distribution is identical to the
desired target distribution:

.. autoclass:: epifx.select.TargetAny

The :class:`~TargetPeakMVN` distribution is a multivariate normal distribution
for the peak timing and size, as defined by previously-observed peaks:

.. autoclass:: epifx.select.TargetPeakMVN

Observation models
------------------

The ``epifx.obs`` module provides generic
:doc:`observation models <observation_models>` for count data with known or
unknown denominators, as well as functions for reading observations from disk.

Each observation model must have a unique ``unit``, and is used to calculate
likelihoods for all observations that share this same ``unit``.

.. code-block:: python

   import epifx.obs

   # Create the simulation parameters.
   params = ...
   # Create an observation model for weekly data (with a period of 7 days),
   # that pertains to all observations whose unit is "obs_unit".
   obs_model = epifx.obs.PopnCounts("obs_unit", obs_period=7)
   # Define the observation model parameters.
   obs_model.define_params(params, bg_obs=300, pr_obs=0.01, disp=100)


Forecast summaries
------------------

.. autofunction:: epifx.summary.make

.. autoclass:: epifx.summary.PrOutbreak

.. autoclass:: epifx.summary.PeakMonitor
   :no-members:
   :members: peak_size, peak_time, peak_date, peak_weight, expected_obs, days_to

.. autoclass:: epifx.summary.PeakForecastEnsembles

.. autoclass:: epifx.summary.PeakForecastCIs

.. autoclass:: epifx.summary.PeakSizeAccuracy

.. autoclass:: epifx.summary.PeakTimeAccuracy

.. autoclass:: epifx.summary.ExpectedObs

.. autoclass:: epifx.summary.ObsLikelihood

.. autoclass:: epifx.summary.ThresholdMonitor
   :no-members:
   :members: exceed_date, exceed_mask, exceed_weight

.. autoclass:: epifx.summary.ExceedThreshold

Generating forecasts
--------------------

.. code-block:: python

   import datetime
   import epifx
   import epifx.obs

   # Simulation parameters
   num_px = 3600
   model = epifx.SEIR()
   time = pypfilt.Datetime()
   popn = 4000000
   prng_seed = 42
   params = epifx.default_params(num_px, model, time, popn_size, prng_seed)

   # Simulate from the 1st of May to the 31st of October, 2015.
   start = datetime.datetime(2015, 5, 1)
   until = start + datetime.timedelta(days=183)

   # Load the relevant observations.
   obs_list = ...
   # Create an observation model for weekly data (with a period of 7 days).
   obs_model = epifx.obs.PopnCounts("obs_unit", obs_period=7)
   # Define the observation model parameters.
   obs_model.define_params(params, bg_obs=300, pr_obs=0.01, disp=100)

   # Generate weekly forecasts for the first 9 weeks (and the start date).
   fs_dates = [start + datetime.timedelta(days=week * 7)
               for week in range(10)]

   # Summary statistics and output file.
   summary = epifx.summary.make(params, obs_list))
   out = "forecasts.hdf5"

   # Generate the forecasts and save them to disk.
   pypfilt.forecast(params, start, until, [obs_list], fs_dates, summary, out)

Comparing forecast outputs
--------------------------

Output files can be compared for equality, which is useful for ensuring that
different systems produce identical results.

.. autofunction:: epifx.cmd.cmp.files

This functionality is also provided as a :ref:`command-line script <cmp>`.
