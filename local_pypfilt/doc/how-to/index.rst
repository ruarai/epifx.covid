.. _how_to:

How-to Guides
=============

**TODO:** provide worked examples of specific tasks.

Creating your own model
-----------------------

Creating your own observation model
-----------------------------------

Defining prior distributions
----------------------------

Provide models with lookup tables
---------------------------------

Provide observation models with lookup tables
---------------------------------------------

An observation model **may** allow some of its parameters to be defined in
lookup tables, rather than as fixed values.
This allows parameters to differ between particles and to vary over time.

.. note:: The observation model must support using lookup tables for parameter
   values.
   The example shown here uses the ``epifx.obs.PopnCounts`` observation model,
   which allows the observation probability to be defined in a lookup table.

To use a lookup table for an observation model parameter, the scenario must:

* Define the lookup table by giving it a name (in this example, ``"pr_obs"``)
  and identifying the data file to use (in this example, ``"pr-obs.ssv"``);
* Notify pypfilt_ that each particle should be associated with a column from
  this lookup table, by adding the table's name to the
  ``"sample_lookup_tables"`` list; and
* Notify the observation model that it should use this lookup table for the
  observation probability, by providing it as an argument to the observation
  model's constructor (in this example, the argument name is
  ``"pr_obs_lookup"``).

These steps are illustrated in the following TOML_ excerpt, which shows only
the relevant lines:

.. code:: toml

   [scenario.test]
   sample_lookup_tables = ["pr_obs"]

   [scenario.test.lookup_tables]
   pr_obs = "pr-obs.ssv"

   [scenario.test.observations.cases]
   model = "epifx.obs.PopnCounts"
   init.pr_obs_lookup = "pr_obs"

Using lookup tables in your own models
--------------------------------------

Using lookup tables in your own observation models
--------------------------------------------------

Defining and using parameters
-----------------------------

Particle filter: resampling
---------------------------

Particle filter: post-regularisation
------------------------------------

Creating summary tables
-----------------------

Creating summary monitors
-------------------------
