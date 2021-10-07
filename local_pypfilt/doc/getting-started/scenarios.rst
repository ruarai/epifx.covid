Forecast scenarios
==================

The example file ``predation.toml`` (shown below) defines all of details that
are needed to generate forecasts.
See the :ref:`concepts` for further information.

.. note:: A single TOML_ file may contain multiple forecast scenarios.
   These can be used to generate separate forecasts for different time
   periods, locations, observation models, etc.

#. The key components are defined in the ``[components]`` table:

   * The simulation model (``model``);
   * The time scale (``time``, either scalar time or date-time); and
   * The output recorder (``summary``).

#. Prior distributions for each model parameter are defined in the
   ``[model.priors]`` table.

#. Lower and upper bounds for each model parameter are defined in the
   ``[model.bounds]`` table.

#. Summary statistics to be recorded in the output file are defined in the
   ``[summary.monitors]`` and ``[summary.tables]`` tables.

#. Various particle filter settings are defined in the ``[parameters]`` table,
   including:

   * The number of particles (``particles``);
   * The seed for the pseudo-random number generator (``prng_seed``); and
   * The simulation time period (``time.start`` and ``time.end``).

It also defines a single **forecast scenario** in the ``[scenario.example]``
table, which contains:

* A unique name to identify this scenario (``"Example Scenario"``);
* An observation model for :math:`x(t)` observations
  (``[scenario.example.observations.x]`` table); and
* The observation model for :math:`y(t)` observations
  (``[scenario.example.observations.y]`` table).

.. literalinclude:: ../../src/pypfilt/examples/predation.toml
   :language: toml
   :linenos:
   :name: predation-toml
   :caption: The contents of ``predation.toml``.
