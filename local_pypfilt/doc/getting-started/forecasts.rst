.. _gs-forecast:

Running the forecasts
=====================

With the forecast scenarios defined in the file ``predation.toml``, generating
the forecasts is as simple as defining the times at which to generate a
forecast, and looping over each forecast scenario in the TOML_ file:

.. code-block:: python

   import pypfilt

   config_file = 'predation.toml'
   forecast_times = [1.0, 3.0, 5.0, 7.0, 9.0]
   for forecast in pypfilt.forecasts_iter(config_file):
       pypfilt.forecast(forecast.params, forecast.observation_streams,
                        forecast_times, filename='output.hdf5')

This will save all of the summary tables for each forecast in the file
``'output.hdf5'``.

.. note:: `HDF5 <http://hdfgroup.org/>`__ is a file format that allows you to
   store lots of data tables and related metadata in a single file, and to
   load these data tables as if they were NumPy arrays.
   All of the summary tables recorded by pypfilt_ are NumPy
   `structured arrays <https://numpy.org/doc/stable/user/basics.rec.html>`__.
   You can explore HDF5 files with the `h5py <https://www.h5py.org/>`__
   package, which makes it easy to load and store data tables.

We'll now look at how to load and plot these outputs.
