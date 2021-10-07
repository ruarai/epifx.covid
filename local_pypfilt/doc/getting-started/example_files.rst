Example files
=============

The ``pypfilt.examples.predation`` module contains all of the code required to
generate forecasts for this model, and the following example files:

``predation-counts-x.ssv``
   Example observations for the prey population :math:`x(t)`.

``predation-counts-y.ssv``
   Example observations for the predator population :math:`y(t)`.

``predation.toml``
   A TOML_ file that defines everything needed to generate forecasts from
   these example data files.

You can easily make local copies of these files:

.. code:: python

   import pypfilt.examples.predation

   pypfilt.examples.predation.write_example_files()

We'll now look at how to define forecast scenarios.
