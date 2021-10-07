pypfilt.time
============

.. py:module:: pypfilt.time

Two pre-defined simulation time scales are provided.

.. autoclass:: pypfilt.time.Scalar
   :no-members:
   :members: __init__, set_period, with_observations,
      with_observations_from_time

.. autoclass:: pypfilt.time.Datetime
   :no-members:
   :members: __init__, set_period, with_observations,
      with_observations_from_time

Custom time scales
------------------

If neither of the above time scales is suitable, you can define a custom time
scale, which should derive the following base class and define the methods
listed here:

.. autoclass:: pypfilt.time.Time
   :exclude-members: stream_tuple, update_tuple, set_period,
      with_observations, with_observations_from_time
