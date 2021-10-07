0.6.0 (2020-08-12)
------------------

This release introduces major structural changes to the entire package, and
incorporates a number of features that were originally implemented in the
``epifx`` package.
Please see the online documentation for further details.
The major user-facing changes are:

* Breaking change: drop support for Python 2, require Python 3.6 or newer.

* Breaking change: forecast scenarios are now defined in TOML files.

0.5.5 (2019-11-25)
------------------

* Bug fix: ensure that ``pypfilt.step`` records the true start of the
  simulation period, if it has not already been defined.

* Enhancement: ``pypfilt.run`` now returns the current index into the history
  matrix, which allows repeat calls to ``pypfilt.run`` to be chained together.
  This may be of use when, e.g., generating a sequence of forecasts where each
  forecast is sufficiently short that it will not cause the simulation window
  to move past the end of the previous estimation run.

* Ensure the documentation builds correctly on Read The Docs.

0.5.4 (2017-10-26)
------------------

* Bug fix: ensure the true start of the simulation period is always recorded.

0.5.3 (2017-10-26)
------------------

* Enhancement: record the true start of the simulation period, so that even if
  the estimation run or forecasting run begins at a later date, the true start
  is available (``params['epoch']``).

* Enhancement: axis and series labels can now be defined by arbitrary
  functions.

* Enhancement: ``pypfilt.plot.series`` now support string scales.

* Enhancement: the ``pypfilt.check`` module provides convenience functions for
  checking invariants. Currently, it is able to check the history matrix
  dimensions. See the API documentation for further details.

* Enhancement: add instructions for install ``pypfilt`` with pip.

* Enhancement: provide example commands for the release process.

0.5.2 (2017-05-05)
------------------

* Bug fix: make ``pypfilt.examples`` a valid Python module.

* Bug fix: fix the Lotka-Volterra model in ``pypfilt.examples.predation`` to
  work correctly with scalar and non-scalar time scales.

0.5.1 (2017-04-28)
------------------

* Bug fix: correctly generate summaries for the case where no table rows will
  be generated. This bug was introduced in pypfilt 0.5.0 (commit ``8a0a614``).

0.5.0 (2017-04-26)
------------------

* Breaking change: the base model class has been renamed to ``pypfilt.Model``.

* Breaking change: the base model class has been simplified; the
  ``state_info``, ``param_info``, and ``param_bounds`` methods have been
  replaced by a single method, ``describe``. This method *also* defines, for
  each element of the state vector, whether that element can be sampled
  continuously (e.g., by the post-regularised filter).

* Breaking change: ``pypfilt.summary.HDF5`` no longer creates a table of
  observations if no such table has been defined, since it may be desirable to
  store observations in multiple tables (e.g., grouped by source or
  observation unit). To retain the previous behaviour, add the new
  observations table ``pypfilt.summary.Obs`` to the summary object.

* Breaking change: particle weights are now passed as an additional argument
  to the log-likelihood function. Previously, the log-likelihood function was
  inspected to determine whether it accepted an extra argument (a nasty hack).

* Bug fix: avoid raising an exception when ``regularise_or_fail`` is ``False``
  (this was the intended behaviour in previous versions).

* Bug fix: ensure that ``pypfilt.summary.obs_table`` correctly encodes the
  observation source and units.

* Bug fix: correct an off-by-one error in ``pypfilt.stats.qtl_wt`` that caused
  the weighted quantiles to be calculated incorrectly. The calculation error
  was inversely proportional to the number of particles and should be
  negligible for any reasonable number of particles (e.g., :math:`\ge 10^3`).

* Enhancement: custom simulation time scales are supported. Two time scales
  are provided (``pypfilt.Datetime`` and ``pypfilt.Scalar``) and additional
  time scales can be implemented by inheriting from ``pypfilt.time.Time``.

* Enhancement: allow likelihoods to depend on past states by settings
  ``params['last_n_periods']`` to N > 1, so that the current observation
  period can be compared to previous observation periods.

* Enhancement: monitor states are now cached and restored, allowing them to
  calculate statistics over the *combined* estimation *and* forecasting runs.
  This means that, e.g., peak times and sizes are correctly reported *even* if
  they occurred prior to the forecasting date.

* Enhancement: add conversion functions for manipulating individual columns in
  structured arrays.

* Enhancement: plotting functions are provided by a new module,
  ``pypfilt.plot`` (adding an *optional* dependency on matplotlib).

* Enhancement: provide a base class for simulation metadata
  (``pypfilt.summary.Metadata``).

* Enhancement: the (continuous) Lotka-Volterra equations are provided as an
  example in ``pypfilt.examples.predation`` and act as the example system in
  the documentation.

* Enhancement: ``pypfilt.summary.dtype_names_to_str`` now also accepts fields
  as a list field names (i.e., strings).

* Enhancement: test cases for several modules are now provided in ``./tests``
  and can be run with `tox <https://tox.readthedocs.io/>`__.

* Enhancement: document how to install required packages as wheels, avoiding
  lengthy compilation times.

* Enhancement: document the release process and provide instructions for
  uploading packages to PyPI.

0.4.3 (2016-09-16)
------------------

* Bug fix: correct the basic resampling method. Previously, random samples
  were drawn from the unit interval and were erroneously assumed to be in
  sorted order (as is the case for the stratified and deterministic methods).

* Enhancement: automatically convert Unicode field names to native strings
  when using Python 2, to prevent NumPy from throwing a
  `TypeError <https://github.com/numpy/numpy/issues/2407>`__, as may occur
  when using ``from __future__ import unicode_literals``.

  This functionality is provided by ``pypfilt.summary.dtype_names_to_str``.

* Enhancement: ensure that temporary files are deleted when the simulation
  process is terminated by the SIGTERM signal.

  Previously, they were only deleted upon normal termination (as noted in the
  `atexit <https://docs.python.org/2/library/atexit.html>`__ documentation).

* Enhancement: consistently separate Unicode strings from bytes, and provide
  utility functions in the ``pypfilt.text`` module.

* Enhancement: forecast from the most recent known-good cached state, avoiding
  the estimation pass whenever possible.

* Enhancement: allow the observation table to be generated externally. This
  means that users can include additional columns as needed.

* Enhancement: separate the calculation of log-likelihoods from the adjustment
  of particle weights, resulting in the new function ``pypfilt.log_llhd_of``.

* Enhancement: provide particle weights to the log-likelihood function, if the
  log-likelihood function accepts an extra argument. This has no impact on
  existing log-likelihood functions.

* Enhancement: by default, allow simulations to continue if regularisation
  fails. This behaviour can be changed::

      params['resample']['regularise_or_fail'] = True

0.4.2 (2016-06-16)
------------------

* Breaking change: ``pypfilt.forecast`` will raise an exception if no
  forecasting dates are provided.

* Add installation instructions for Red Hat Enterprise Linux, Fedora, and Mac
  OS X (using `Homebrew <http://brew.sh/>`__).

0.4.1 (2016-04-26)
------------------

* Enhancement: allow forecasts to resume from cached states, greatly improving
  the speed with which forecasts can be generated when new or updated
  observations become available. This is enabled by defining a cache file::

      params['hist']['cache_file'] = 'cache.hdf5'

* Enhancement: add option to restrict summary statistics to forecasting
  simulations, ignoring the initial estimation run. This is enabled by passing
  ``only_fs=True`` as an argument to the ``pypfilt.summary.HDF5`` constructor.

0.4.0 (2016-04-22)
------------------

* Breaking change: require models to define default parameter bounds by
  implementing the ``param_bounds`` method.

* Enhancement: offer the post-regularised particle filter (post-RPF) as an
  alternative means of avoiding particle impoverishment (as opposed to
  incorporating stochastic noise into the model equations). This is enabled by
  setting::

      params['resample']['regularisation'] = True

  See the example script (``./doc/example/run.py``) for a demonstration.

* Improved documentation for ``pypfilt.model.Base`` and summary statistics.

* Add documentation for installing in a virtual environment.

0.3.0 (2016-02-23)
------------------

* This release includes a complete overhaul of simulation metadata and summary
  statistics. See ``./doc/example/run.py`` for an overview of these changes.

* Breaking change: decrease the default resampling threshold from 75% to 25%.

* Breaking change: define base classes for summary statistics and output.

* Breaking change: define a base class for simulation models.

* Breaking change: collate the resampling and history matrix parameters to
  reduce clutter.

* Breaking change: move ``pypfilt.metadata_priors`` to ``pypfilt.summary``.

* Bug fix: prevent ``stats.cov_wt`` from mutating the history matrix.

* Bug fix: ensure that the time-step mapping behaves as documented.

* Bug fix: ensure that state vector slices have correct dimensions.

* Enhancement: ensure that forecasting dates lie within the simulation period.

* Performance improvement: Vectorise the history matrix initialisation.

* Host the documentation at Read The Docs.

0.2.0 (2015-11-16)
------------------

* Notify models whether the current simulation is a forecast (i.e., if there
  are no observations). This allows deterministic models to add noise when
  estimating, to allow identical particles to differ in their behaviour, and
  to avoid doing so when forecasting.

  Note that this is a breaking change, as it alters the parameters passed to
  the model update function.

* Simplify the API for running a single simulation; ``pypfilt.set_limits`` has
  been removed and ``pypfilt.Time`` is not included in the API documentation,
  on the grounds that users should not need to make use of this class.

* Greater use of NumPy array functions, removing the dependency on six >= 1.7.

* Minor corrections to the example script (``./doc/example/run.py``).

0.1.2 (2015-06-08)
------------------

* Avoid error messages if no logging handler is configured by the application.

* Use a relative path for the output directory. This makes simulation metadata
  easier to reproduce, since the absolute path of the output directory is no
  longer included in the output file.

* Build a universal wheel via ``python setup.py bdist_wheel``, which supports
  both Python 2 and Python 3.


0.1.1 (2015-06-01)
------------------

* Make the output directory a simulation parameter (``out_dir``) so that it
  can be changed without affecting the working directory, and vice versa.


0.1.0 (2015-05-29)
------------------

* Initial release.
