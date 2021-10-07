0.6.1 (2020-10-20)
------------------

* Bug fix: handle peak size and time statistics for peaks of zero cases.

0.6.0 (2020-08-13)
------------------

This release introduces major structural changes to the entire package. Please
see the online documentation for further details. The major user-facing
changes are:

* Breaking change: drop support for Python 2, require Python 3.6 or newer.

* Breaking change: forecast scenarios are now defined in TOML files.

* Breaking change: use the new ``epifx-decl-fs`` command to run forecasts.

* Breaking change: the previous ``epifx-xxx`` commands have not been updated
  to work with the new TOML forecast scenarios, and are not currently
  supported. They will be updated in a future release.

0.5.8 (2020-04-03)
------------------

* Enhancement: ``epifx.obs.PopnCount`` now supports incomplete observations
  that comprise an incomplete count and an estimated upper bound for the true
  value.

0.5.7 (2020-03-28)
------------------

* Enhancement: add an SEEIIR model, which has the same parameters as the SEIR
  model.

0.5.6 (2018-06-07)
------------------

* Enhancement: observation upper bounds can now be treated as point estimates
  by setting ``params['epifx']['upper_bounds_as_obs'] = True``.

* Enhancement: ``epifx-replay`` can now define "perfect" upper bound estimates
  (corresponding to the observed values in the most recent data snapshot) by
  using the ``--perfect-upper-bounds`` argument.

* Enhancement: ``epifx.cmd.run_in_parallel`` now returns the job completion
  status (Boolean value).

* Test cases are now run against Python 3.6 rather than Python 3.5, since
  Debian Testing has migrated to Python 3.6.

0.5.5 (2017-10-26)
------------------

* This release requires pypfilt >= 0.5.4.

* Bug fix: correctly calculate the time of first infection relative to the
  true start of the simulation period (``params['epoch']``) rather than, e.g.,
  the start of the forecasting run.

* Bug fix: ignore duplicate table rows when generating JSON files.

* Bug fix: additional data validation checks when generating JSON files.

* Enhancement: add a new command (``epifx-replay``) that replays the
  observations from an existing set of forecasts, so that new forecasts can be
  generated while accounting for, e.g., incomplete observations.

* Enhancement: add a worked example of generating forecasts. See the "Worked
  Example" page in the online documentation.

0.5.4 (2017-08-17)
------------------

* Bug fix: when processing forecast files that only contain estimation runs,
  the ``epifx-json`` command will ignore forecast credible intervals prior to
  the date of the most recent observation.

0.5.3 (2017-08-16)
------------------

* Bug fix: ensure ``epifx-json`` uses native strings for dtype field names.

* Bug fix: remove invalid import statements that rely on an as-yet unreleased
  version of ``pypfilt``.

0.5.2 (2017-08-16)
------------------

* Enhancement: add a new monitor, ``epifx.summary.ThresholdMonitor``, and a
  corresponding table, ``epifx.summary.ExceedThreshold``, for detecting when
  the expected observations for each particle exceed a specific threshold.
  Note that this table is **not** included in the default summary object
  returned by ``epifx.summary.make``.

* Enhancement: the ``epifx-json`` command no longer aborts if a forecast file
  only contains forecast credible intervals (the ``forecasts`` table) but not
  the peak timing credible intervals (the ``peak_cints`` table).

* Enhancement: the ``epifx-json`` command will generate output from estimation
  runs if a forecast file only contains an estimation run (which can occur,
  e.g., when directly using ``pypfilt.run`` to generate forecasts rather than
  using ``pypfilt.forecast``).

0.5.1 (2017-06-15)
------------------

* Bug fix: the accept-reject sampler now resets particle weights after each
  iteration. This only affects summary tables that require the weights to sum
  to unity, it has no effect on the particle selection.

* Bug fix: the daily forcing signal now correctly returns datetime instances.

* Bug fix: add a missing function argument when processing incomplete
  observations.

* Enhancement: the ``epifx-forecast`` command now accepts a new argument,
  ``--all`` (or ``-a``), which generates forecasts for all defined locations.

0.5.0 (2017-05-10)
------------------

* This release requires pypfilt >= 0.5.1.

* Breaking change: ``epifx.default_params`` makes fewer assumptions, and now
  takes more positional arguments.

* Breaking change: the SEIR model is now defined in terms of intrinsic
  parameters (e.g., R0 rather than the daily force of infection) and the time
  of initial exposure is now just another model parameter.

* Breaking change: record predictive CIs in the ``/data/forecasts`` table; by
  default, the median and the 50% and 95% credible intervals are recorded. The
  expected observation CIs are now stored in the ``/data/expected_obs`` table.

* Breaking change: record observations in tables grouped by the observation
  unit; these tables are now located in HDF5 tables ``/data/obs/obs_unit``.

* Bug fix: ensure that ``epifx.summary.ObsLikelihood`` correctly encodes the
  observation source and units.

* Bug fix: robustly handle near-zero probability mass.

* Enhancement: arbitrary model priors are supported via the accept-reject
  sampler provided by the ``epifx.select`` module.

* Enhancement: a suite of commands for performing retrospective observation
  model scans and live forecasts are provided by the ``epifx.cmd`` module. See
  the documentation for details.

* Enhancement: an example template is provided (see the ``epifx.example``
  module and the ``epifx-template`` command) that includes Australian Google
  Flu Trends data. Observations can be loaded with ``epifx.example.gft_obs``.

* Enhancement: ``epifx.summary.PeakMonitor`` now surveys the entire particle
  trajectory (including both the estimation and forecasting runs) so that
  peaks that occurred prior to the forecasting date are reported correctly.

* Enhancement: the ``epifx.summary.ObsLikelihood`` table can now record the
  likelihood of arbitrary observations (i.e., observations that were **not**
  included in the filtering process).

* Enhancement: the default summary tables provided by ``epifx.summary.make``
  can be suppressed as needed.

* Enhancement: custom simulation time scales are supported.

* Enhancement: add quantile and probability mass sum functions to the
  observation models.

* Enhancement: test cases for several modules are now provided in ``./tests``
  and can be run with `tox <https://tox.readthedocs.io/>`__.

* Enhancement: document the release process and provide instructions for
  uploading packages to PyPI.

0.4.3 (2016-09-16)
------------------

* This release requires pypfilt >= 0.4.3.

* Breaking change: the ``epifx.obs.SampleCounts`` observation model now uses a
  *Beta-binomial* distribution rather than a *Beta* distribution. Parameter
  names and definitions have been changed accordingly.

* Enhancement: consistently separate Unicode strings from bytes, and
  automatically convert NumPy field names into native strings.

* Enhancement: add support for incomplete data for which there may or may not
  be an upper bound (whether known in fact or estimated).

* Enhancement: record the likelihood of each observation according to each
  particle (see the ``epifx.summary.ObsLikelihood`` class).

0.4.2 (2016-06-16)
------------------

* Breaking change: replace the observation models added in epifx 0.4.1 with
  observations models for:

  * Count data where the denominator is *known or assumed* to be the entire
    population (``epifx.obs.PopnCounts``); and

  * Count data where the denominator is reported and may vary, and where the
    background signal is a *fixed proportion* (``epifx.obs.SampleCounts``).

0.4.1 (2016-04-22)
------------------

* Enhancement: provide generic negative binomial observation models for count
  data and for fractional/percentage data in ``epifx.obs``.

0.4.0 (2016-04-22)
------------------

* This release requires pypfilt >= 0.4.0.

* Breaking change: models must define default parameter bounds by implementing
  the ``param_bounds`` method.

* Breaking change: model expectation functions now receive the previous and
  current state vectors, in addition to the infection probability vector. This
  means that expectation functions will need to change from::

      expect(params, unit, period, pr_inf)

  to::

      expect(params, unit, period, pr_inf, prev, curr)

* Enhancement: ``epifx.summary.make`` now passes additional keyword arguments
  to the ``pypfilt.summary.HDF5`` constructor, allowing users to override
  default settings, such as ``first_day=False``.

* Bug fix: ensure that infection probabilities are strictly non-negative.

* Bug fix: ensure that population invariants are enforced correctly.

* Bug fix: correctly scale the daily seeding probability.

* Add instructions for installing epifx in a virtual environment.

0.3.1 (2016-02-25)
------------------

* Bug fix: prevent a runtime error with ``params['epifx']['one_prng'] = True``
  by correctly retrieving the pypfilt PRNG (``params['resample']['rnd']``).

0.3.0 (2016-02-23)
------------------

* This release requires pypfilt >= 0.3.0.

* Provide each summary statistic as a separate class.

* Inherit from the pypfilt simulation model base class.

* Host the documentation at Read The Docs.

0.2.0 (2015-11-17)
------------------

* Use an independent PRNG instance for model stochasticity, as distinct from
  the pypfilt PRNG instance (used for resampling). Note that this breaks
  backwards compatibility (in terms of producing identical outputs) and can be
  overridden by setting ``params['epifx']['one_prng'] = True``.

* This release requires pypfilt >= 0.2.0.

0.1.10 (2015-07-09)
-------------------

* Fix a bug where temporal forcing could cause a negative force of infection.


0.1.9 (2015-07-06)
------------------

* Add support for temporal forcing, modulated by the new parameter sigma.


0.1.8 (2015-06-18)
------------------

* Update the model parameter invariants for alpha, based on the priors for R0
  and gamma.


0.1.7 (2015-06-18)
------------------

* Sample R0 and calculate alpha, rather than sampling alpha directly.


0.1.6 (2015-06-08)
------------------

* Avoid error messages if no logging handler is configured by the application.

* Default to comparing only the simulation outputs and ignore the metadata;
  this can be overridden by the ``--meta-data`` (``-m``) option.

* Build a universal wheel via ``python setup.py bdist_wheel``, which supports
  both Python 2 and Python 3.

* This release requires pypfilt >= 0.1.2.


0.1.5 (2015-06-04)
------------------

* Record credible intervals for state variables (``/data/state``).


0.1.4 (2015-06-03)
------------------

* Reduce the minimum latent period to half a day.

* No longer require the simulation period to be included in the simulation
  parameters dictionary.

* Hide stderr output from spawned processes when obtaining git metadata, since
  the error messages have no relevance to the user (they only serve to
  indicate that the working directory is not part of a git repository).


0.1.3 (2015-06-01)
------------------

* Obtain git metadata from the working directory, if it is contained within a
  repository. This requires pypfilt >= 0.1.1.

  Note that sensible metadata will only be obtained if the working directory
  is *not* manipulated by program code (e.g., ``os.chdir``).


0.1.2 (2015-06-01)
------------------

* Record the enforced limits on model parameters, so that output files include
  sufficient information to independently produce identical results.

  The default limits on ``beta`` and ``gamma`` are now identical to the domain
  of their default priors (1 to 3 days).

* Added options ``--verbose`` and ``--data-only`` to the ``cmp-output``
  script.

* Ignore the command line (``/meta/sim/cmdline``) when comparing output files.


0.1.1 (2015-05-29)
------------------

* Added a script (``cmp-output``) that compares output files for identical
  simulation outputs and metadata.


0.1.0 (2015-05-29)
------------------

* Initial release.
