pypfilt.io
==========

.. py:module:: pypfilt.io

The :mod:`pypfilt.io` module provides functions for reading tabular data from
data files.

Reading tabular data
--------------------

.. autofunction:: pypfilt.io.read_table

.. autofunction:: pypfilt.io.date_column

.. autofunction:: pypfilt.io.datetime_column

Lookup tables
-------------

The :mod:`pypfilt.io` module also provides lookup tables, which are used to
retrieve time-indexed values (e.g., time-varying model inputs).

.. autofunction:: pypfilt.io.read_lookup_table

.. autofunction:: pypfilt.io.lookup_values_count

.. autofunction:: pypfilt.io.lookup
