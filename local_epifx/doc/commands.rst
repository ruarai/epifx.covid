Commands
========

The following commands are part of `epifx`_:

- :literal:`epifx-template` creates an example :ref:`configuration <config>`,
  which serves as a starting point for defining data sources, observation
  models, and simulation settings.

- :literal:`epifx-locns` lists the locations that are defined in the
  :ref:`local configuration <config>`.

- :literal:`epifx-scan` generates forecasts from retrospective data, over a
  range of parameter values for the observation model(s), so that optimal
  observation model settings can be :ref:`identified <scan>`.

- :literal:`epifx-summary` :ref:`summarises <summary>` forecasting performance
  for a completed observation model scan.

- :literal:`epifx-forecast` generates :ref:`forecasts <forecast>` from live
  data (i.e., where the future observations are unknown).

- :literal:`epifx-replay` will :ref:`replay <replay>` the observations from an
  existing set of forecasts, to produce a new set of forecasts (e.g.,
  accounting for incomplete observations).

- :literal:`epifx-json` converts forecasts into JSON files that can be plotted
  and :ref:`viewed interactively <json>` in a web browser.

- :literal:`epifx-cmp` compares simulation outputs between files, which is
  useful for ensuring that different systems produce
  :ref:`identical results <cmp>`.

.. _config:

Local configuration
-------------------

All of the `epifx`_ commands rely on the necessary simulation settings being
defined in a user-defined module, ``epifxlocns``.
Run the follow command to generate a basic template for this module in the
current directory:

.. code-block:: shell

   epifx-template

The :literal:`epifx-locns` command lists the locations that are defined in
this module:

.. code-block:: shell

   epifx-locns

This module must define the following function:

.. function:: local_settings(locn_id=None)

   Return the simulation settings for the user-defined location uniquely
   identified by ``locn_id`` or, if ``locn_id is None``, return a list of the
   user-defined locations.

The settings for each location should be returned as a dictionary that
contains the following keys:

- ``'id'``: the unique identifier for this location (**should not** contain
  spaces or numbers).
- ``'name'``: a pretty-printed version of the location identifier (**may**
  contain spaces and numbers).
- ``'popn'``: the population size.
- ``'obs_model'``: an :ref:`observation model <obsmodels>` instance.
- ``'obs_file'``: the name of the observations file.
- ``'obs_filter'``: a function use to remove outliers from the observations,
  if required (otherwise, this key is not required).
- ``'from_file_args'``:  a dictionary of additional keyword arguments to pass
  to the ``from_file`` method of the observation model (e.g.,
  ``'value_col'`` for :func:`epifx.obs.PopnCounts.from_file`).
- ``'scan_years'``: a list of the years for which observation model scans
  should be performed.
- ``'scan'``: a dictionary that maps observation model parameter names to one
  or more values for that parameter, to be covered by observation model scans.
  Values may be defined as scalar values, lists, or dictionaries that map
  years to either scalar values or lists.
- ``'forecast'``: a dictionary that maps observation model parameter names to
  one or more values for that parameter, to be used when forecasting.
  Values may be defined as scalar values, lists, or dictionaries that map
  years to either scalar values or lists.
- ``'om_format'``: a dictionary that maps observation model parameter names to
  format specifiers, which are used to include parameter values in output file
  names.
- ``'om_name'``: a dictionary that maps observation model parameter names to
  strings, which are used to identify these parameters in output file names.
- ``'out_dir'``: the directory in which simulation outputs will be written.
- ``'json_dir'``: the directory in which JSON forecast files will be
  written.
- ``'tmp_dir'``: the directory in which temporary files will be stored when
  performing simulations.
- ``'get_params'``: a function that accepts this dictionary as an argument and
  returns the simulation parameters dictionary.
- JSON plot settings:

  - ``'obs_axis_lbl'``: the axis label for the observations.
  - ``'obs_axis_prec'``: the decimal precision of axis ticks for the
    observations.
  - ``'obs_datum_lbl'``: the label for individual observations.
  - ``'obs_datum_prec'``: the decimal precision of individual observations.

- ``'extra_args'``: used to apply custom settings arguments for observation
  model scans and for forecasts:

  - ``'start'``: a function that takes one argument, the season, and returns
    the start of the simulation period.
  - ``'until'``: a function that takes one argument, the season, and returns
    the end of the simulation period.
  - ``'live_fs_dates'``: a function that takes three arguments: the season,
    a list of all observations, and an (optional) initial forecasting date,
    and returns a list of dates for which live forecasts should be generated.
  - ``'scan_fs_dates'``: a function that takes two arguments, the season and
    the list of observations for that season, and returns the dates for which
    retrospective forecasts should be generated.
  - ``'make_summary'``: a function that takes the same arguments as
    :py:func:`epifx.summary.make` and returns a summary object.

.. code-block:: shell

   epifx-template

.. _scan:

Observation model scans
-----------------------

To run simulations for every combination of observation model parameters (as
defined in the :ref:`local configuration <config>`) use :literal:`epifx-scan`.
In the example below, the simulations will be run in 16 parallel processes and
will only be run against observations for the 2015 calendar year:

.. code-block:: shell

   epifx-scan --spawn 16 --year 2015 location_name

.. _summary:

Observation model performance
-----------------------------

To summarise the forecasting performance for an observation model scan, use
:literal:`epifx-summary`:

.. code-block:: shell

   epifx-summary location_name

By default, this will also convert the best forecast for each calendar year
into a JSON file, for interactive viewing in a web browser.

.. _forecast:

Live forecasting
----------------

To generate forecasts from live data, use :literal:`epifx-forecast`.
In the example below, the most recent observation is deemed to be "incomplete"
(i.e., an underestimate of the real value) and an upper bound on the actual
value is provided:

.. code-block:: shell

   epifx-forecast --incomplete 1 --upper-bound 200 location_name

By default, each forecast will also be converted into a JSON file, for
interactive viewing in a web browser.

.. _replay:

Replay forecasts
----------------

To generate forecasts against existing data, where a given observation may
change over time (e.g., gradual accumulation of count data), use
:literal:`epifx-replay` to run forecasts against the data **as they were
provided** at each forecasting date.
The observations and forecasting dates will be read from an existing JSON
file, as produced by :ref:`epifx-forecast <forecast>` or
:ref:`epifx-json <json>`:

.. code-block:: shell

   epifx-replay prev-forecasts.json location_name

By default, the forecasting location specified in the JSON file will also be
used to generate the new forecasts, but a different forecasting location can
be specified (as illustrated above).

Incomplete observations can also be adjusted to have "perfect" upper bound
estimates (i.e., corresponding to the observed values in the most recent data
snapshot):

.. code-block:: shell

   epifx-replay --perfect-upper-bounds prev-forecasts.json location_name

.. note::

   If the final reported value for an observation date is **smaller** than the
   reported value in an earlier snapshot, the original observation will be
   left **unchanged**.

   This may be changed in a future version of epifx_, but will require
   changing the semantics (and probably the name) of the observation
   ``'upper_bound'`` field and the observation model log-likelihood functions.

.. _json:

Interactive forecast plots
--------------------------

Both :literal:`epifx-summary` and :literal:`epifx-forecast` will produce JSON
files by default, but :literal:`epifx-json` can be used to convert **any** set
of forecast output files into a JSON file:

.. code-block:: shell

   epifx-json --location location_name --output forecasts.json *.hdf5

.. _cmp:

Comparing simulation outputs
----------------------------

Output files can be compared for equality, which is useful for ensuring that
different systems produce identical results.

.. code-block:: shell

   epifx-cmp --help
   epifx-cmp file1.hdf5 file2.hdf5

Note that simulation outputs have been observed to differ depending on the
installed versions of the NumPy and SciPy packages, due to changes in the
random number generator.
