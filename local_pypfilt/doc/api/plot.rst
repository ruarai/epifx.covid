.. _api_plot:

pypfilt.plot
============

.. py:module:: pypfilt.plot

Several plotting routines, built on top of
`matplotlib <http://matplotlib.org/>`__, are provided in the
:mod:`pypilt.plot` module.

.. note:: `matplotlib <http://matplotlib.org/>`__ **must** be installed in
          order to use this module.

To generate plots non-interactively (i.e., without having a window appear) use
the ``'Agg'`` backend::

    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

See the
`matplotlib FAQ <http://matplotlib.org/faq/howto_faq.html#generate-images-without-having-a-window-appear>`__ for more details.

Styles and colour palettes
--------------------------

.. autofunction:: default_style

.. autofunction:: apply_style

.. autofunction:: n_colours

.. autofunction:: brewer_qual

.. autofunction:: colour_iter

Plotting functions
------------------

.. autofunction:: cred_ints

.. autofunction:: observations

.. autofunction:: series

Faceted plots
-------------

This package provides a base class (:class:`~Plot`) for plots that comprise
any number of subplots, and three subclasses for specific types of plots:

* :class:`~Wrap` for plots where a single variable identifies each subplot.
* :class:`~Grid` for plots where two variables are used to identify each
  subplot.
* :class:`~Single` for single plots.

The key method of these classes is :func:`Plot.subplots`, which returns an
iterator that yields ``(axes, data)`` tuples for each subplot.
By looping over these tuples, one set of plotting commands can be used to
generate all of the subplots.
For examples, see the ``plot_forecasts()`` and ``plot_params()`` functions in
the :ref:`gs-plotting` section of :ref:`getting_started`.

.. autoclass:: Plot

.. autoclass:: Wrap

.. autoclass:: Grid

For consistency, a class is also provided for single plots.

.. autoclass:: Single
