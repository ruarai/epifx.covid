pypfilt.summary
===============

.. py:module:: pypfilt.summary

Simulation metadata
-------------------

Every simulation data file should include metadata that documents the
simulation parameters and working environment.
The :class:`~pypfilt.summary.Metadata` class provides the means for generating
such metadata:

.. autoclass:: pypfilt.summary.Metadata

Summary data files
------------------

The :class:`HDF5` class encapsulates the process of calculating and recording
summary statistics for each simulation.

.. autoclass:: pypfilt.summary.HDF5
   :no-members:
   :members: save_forecasts

Summary statistic tables
------------------------

Summary statistics are stored in tables, each of which comprises a set of
named columns and a specific number of rows.

The Table class
^^^^^^^^^^^^^^^

To calculate a summary statistic, you need to define a subclass of the
:class:`Table` class and provide implementations of each method.

.. autoclass:: pypfilt.Table
   :members:

Predefined statistics
^^^^^^^^^^^^^^^^^^^^^

The following derived classes are provided to calculate basic summary
statistics of any generic simulation model.

.. autoclass:: pypfilt.summary.ModelCIs

.. autoclass:: pypfilt.summary.ParamCovar
   :special-members:

.. autoclass:: pypfilt.summary.Obs

.. autoclass:: pypfilt.summary.SimulatedObs

.. autoclass:: pypfilt.summary.PredictiveCIs

Utility functions
^^^^^^^^^^^^^^^^^

The following column types are provided for convenience when defining custom
:class:`Table` subclasses.

.. autofunction:: pypfilt.summary.dtype_value

The following functions are provided for converting column types in structured
arrays.

.. autofunction:: pypfilt.summary.convert_cols

.. autofunction:: pypfilt.summary.default_converters

Retrospective statistics
------------------------

In some cases, the :class:`Table` model is not sufficiently flexible, since it
assumes that statistics can be calculated during the course of a simulation.
For some statistics, it may be necessary to observe the entire simulation
before the statistics can be calculated.

In this case, you need to define a subclass of the :class:`Monitor` class,
which will observe ("monitor") each simulation and, upon completion of each
simulation, can calculate the necessary summary statistics.

Note that a :class:`Table` subclass is **also required** to define the table
columns, the number of rows, and to record each row at the end of the
simulation.

.. autoclass:: pypfilt.Monitor
   :members:

Predefined monitors
^^^^^^^^^^^^^^^^^^^

The :class:`PredictiveCIs` summary table requires the following monitor:

.. autoclass:: pypfilt.summary.ExpectedObsMonitor
   :members: expected_obs

Tables and Monitors
-------------------

The methods of each :class:`Table` and :class:`Monitor` will be called in the
following sequence by the :class:`HDF5` summary class:

#. Before any simulations are performed:

   * :func:`Table.dtype`
   * :func:`Monitor.prepare`

   In addition to defining the column types for each :class:`Table`, this
   allows objects to store the simulation parameters and observations.

#. At the start of each simulation:

   * :func:`Monitor.begin_sim`
   * :func:`Table.n_rows`

   This notifies each :class:`Monitor` and each :class:`Table` of the
   simulation period, the number of observation systems (i.e., data sources),
   and whether it is a forecasting simulation (where no resampling will take
   place).

#. During each simulation:

   * :func:`Monitor.monitor`
   * :func:`Table.add_rows`

   This provides a portion of the simulation period for analysis by each
   :class:`Monitor` and each :class:`Table`.
   Because all of the :func:`Monitor.monitor` methods are called before the
   :func:`Table.add_rows` methods, tables can interrogate monitors to obtain
   any quantities of interest that are calculated by :func:`Monitor.monitor`.

#. At the end of each simulation:

   * :func:`Monitor.end_sim`
   * :func:`Table.finished`

   This allows each :class:`Monitor` and each :class:`Table` to perform any
   final calculations once the simulation has completed.
   Because all of the :func:`Monitor.end_sim` methods are called before the
   :func:`Table.finished` methods, tables can interrogate monitors to obtain
   any quantities of interest that are calculated by :func:`Monitor.end_sim`.
