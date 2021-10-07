Epidemic forecasting with mechanistic infection models
======================================================

Description
-----------

This package generates epidemic forecasts with mechanistic infection models.

Documentation
-------------

The documentation is `available online <https://epifx.readthedocs.io/>`_ and
can be built locally with `sphinx <http://sphinx-doc.org/>`_::

    python setup.py build_sphinx

Note that the `sphinx_rtd_theme <https://github.com/snide/sphinx_rtd_theme/>`_
theme must be installed.

License
-------

The code is distributed under the terms of the BSD 3-Clause license (see
``LICENSE``), and the documentation is distributed under the terms of the
`Creative Commons BY-SA 4.0 license
<http://creativecommons.org/licenses/by-sa/4.0/>`_.

Installation
------------

Clone this repository and execute::

    python setup.py install

If you don't have admin rights, install the package locally::

    python setup.py install --user

Dependencies
------------

This package requires `h5py <http://www.h5py.org/>`_ >= 2.2,
`pypfilt <http://bitbucket.org/robmoss/particle-filter-for-python/>`_ >=
0.5.1,
`NumPy <http://www.numpy.org/>`_ >= 1.8, and
`SciPy <http://www.scipy.org/>`_ >= 0.11.
