Testing with tox
================

The pypfilt_ testing suite uses the `pytest <http://docs.pytest.org/>`__
framework, and uses the `tox <https://tox.readthedocs.io/>`__ automation tool
to run the tests under Python 3.
The test cases are contained in the ``./tests`` directory.

To run all tests using all of the Python versions defined in ``tox.ini``, run:

.. code-block:: shell

   tox

The ``tox.ini`` contents are shown below, and include targets that check
whether the documentation in ``./doc`` builds correctly.

.. literalinclude:: ../tox.ini
   :language: ini
