Conclusion
==========

Congratulations!
You've used pypfilt_ to combine a simulation model with observations, generate
forecasts, and plot the results.

You can now try modifying the observation files, adjusting the prior
distributions, etc, and see how this affects the forecast predictions.
Either modify the provided forecast scenario, or add new scenarios to the
example TOML file.
The following script, which reproduces all of the steps in this Getting
Started guide, can be used as a starting point:

.. code:: python

   import pypfilt
   import pypfilt.examples.predation

   pypfilt.examples.predation.write_example_files()

   config_file = 'predation.toml'
   forecast_times = [1.0, 3.0, 5.0, 7.0, 9.0]
   data_file = 'output.hdf5'

   for forecast in pypfilt.forecasts_iter(config_file)
       pypfilt.forecast(forecast.params, forecast.observation_streams,
                        forecast_times, filename=data_file)

   pypfilt.examples.predation.plot(data_file, png=True, pdf=False)

To learn how to use pypfilt_ with your own models and data, see
:ref:`concepts` and the :ref:`how_to`.
