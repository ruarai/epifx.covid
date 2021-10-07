.. _install:

Installation
============

You must have Python 3.6, or later, installed.
On Windows, the simplest option is to use
`Anaconda <https://docs.continuum.io/anaconda/>`__.
By default, it will install all of the required packages except ``toml``.

You can install pypfilt_ with ``pip``. This is best done in a
`virtual environment
<http://docs.python-guide.org/en/latest/dev/virtualenvs/>`__.
It will also install any required packages that are not currently installed.

Install pypfilt without plotting support:
   .. code-block:: shell

      pip install pypfilt

Install pypfilt with plotting support:
   .. code-block:: shell

      pip install pypfilt[plot]

You can also use a package manager, such as ``apt`` (Debian, Ubuntu),
``yum`` (Red Hat Enterprise Linux, CentOS), ``dnf`` (Fedora), or
`Homebrew <http://brew.sh/>`__ (OS X).

The packages required by pypfilt_ are:

* `NumPy <http://www.numpy.org/>`__ 1.17 or newer;
* `SciPy <http://www.scipy.org/>`__ 0.11 or newer;
* `h5py <http://www.h5py.org/>`__ 2.2 or newer;
* `toml <https://github.com/uiri/toml>`__ 0.10 or newer; and
* `matplotlib <http://matplotlib.org/>`__ 1.5 or newer (optional, see
  :ref:`api_plot`).
