Epidemic forecasting in Python
==============================

Welcome to the epifx_ documentation.
This package uses a bootstrap particle filter (pypfilt_) to generate forecasts
of infectious disease epidemics.

.. note:: This documentation assumes that you are already familiar with the
          pypfilt_ package.
          If this is not the case, please read the pypfilt_ **Getting
          Started** tutorial before proceeding.

The epifx_ package provides a number of new components and features on top of
those provided by the pypfilt_ package, including:

* Epidemic models
* :ref:`Observation models <obsmodels>` for case count data
* Epidemic summary statistics
* Fitting flexible model prior distributions
* Running forecasts from the command-line (todo: worked TOML example)

.. _install:

Installation
------------

You can install epifx_ using ``pip``, preferably in a
`virtual environment
<http://docs.python-guide.org/en/latest/dev/virtualenvs/>`__:

   .. code-block:: shell

      pip install epifx

See the `pypfilt documentation <http://pypfilt.readthedocs.io/en/latest/>`__
for further details about installation.

License
-------

The code is distributed under the terms of the BSD 3-Clause license (see
``LICENSE``), and the documentation is distributed under the terms of the
`Creative Commons BY-SA 4.0 license
<http://creativecommons.org/licenses/by-sa/4.0/>`_.

.. _user-docs:

.. toctree::
   :maxdepth: 2
   :caption: User Documentation

   Home <self>
   getting_started
   observation_models
   commands

.. _dev-docs:

.. toctree::
   :maxdepth: 2
   :caption: Development

   contributing
   testing
   release-process
   changelog
