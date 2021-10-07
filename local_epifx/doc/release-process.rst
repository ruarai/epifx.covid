Release process
===============

Feature development takes places on the "master" branch.
Periodically, a release is created by increasing the version number and
tagging the relevant commit with the new version number.

* Check that the release passes all of the tests:

  .. code-block:: shell

     tox

* Update the version number according to the
  `versioning scheme <https://www.python.org/dev/peps/pep-0440/>`__.

   * Update the version number in ``doc/conf.py``.
     The full version must **always** be updated, the short (X.Y) version
     **does not** need to be updated if the version number is being increased
     from X.Y.Z to X.Y.Z+1.

   * Update the version number in ``src/epifx/version.py``.

   * Update the version number in ``setup.cfg``.

* Describe the changes at the top of ``NEWS.rst`` under a heading of the form
  ``X.Y.Z (YYYY-MM-DD)``, which identifies the new version number and the
  date on which this version was released.

* Commit these changes; set the commit message to ``Release epifx X.Y.Z``.

  .. code-block:: shell

     git add NEWS.rst doc/conf.py src/epifx/version.py setup.cfg
     git commit -m "Release epifx X.Y.Z"

* Tag this commit ``X.Y.Z``.

  .. code-block:: shell

     git tag -a X.Y.Z -m "epifx X.Y.Z"

* Push this commit **and** the new tag upstream.

  .. code-block:: shell

     git push --follow-tags

Publishing to PyPI
------------------

These instructions are based on the
`Python Packaging User Guide <https://packaging.python.org/distributing/>`__.

Ensure that ``twine`` is installed:

.. code-block:: shell

   pip install twine

Ensure that all uncommitted changes are stashed, **or they will be packaged!**

.. code-block:: shell

   git stash

Build the wheel ``./dist/epifx-X.Y.Z-py3-none-any.whl``:

.. code-block:: shell

   python setup.py bdist_wheel

Upload this wheel to the PyPI **test** server, so that any problems can be
`identified <https://testpypi.python.org/pypi/epifx/>`__ and fixed:

.. code-block:: shell

   twine upload -r testpypi dist/epifx-X.Y.Z-py3-none-any.whl

Then upload this wheel to PyPI:

.. code-block:: shell

   twine upload dist/epifx-X.Y.Z-py3-none-any.whl

Packaging notes
---------------

To ensure that package data files are included in both source and binary
distributions, (e.g., as created by running ``python setup.by sdist`` and
``python setup.py bdist_wheel``, respectively):

+ Data files should be listed in ``MANIFEST.in``, otherwise they will not be
  included in source distributions.
  Note that these files **will not be included** in binary distributions.

+ By using the setuptools-specific option ``include_package_data`` in
  ``setup.cfg``, when building a binary distribution it will automatically
  include any data files inside the package directories that are specified in
  ``MANIFEST.in``.

  + Note that if you use ``package_data`` in ``setup.cfg`` to specify data
    files, those files **will not be added** unless they are also listed in
    ``MANIFEST.in`` (see the `setuptools`_ documentation).

+ So it appears simplest to define all of the data files in ``MANIFEST.in``
  and set ``include_package_data = True`` in ``setup.cfg``.

.. warning::

   If you update the list of data files, delete the file
   ``src/epifx.egg-info/SOURCES.txt`` before calling ``setup.py`` (see the
   `setuptools`_ documentation).

.. _setuptools: https://setuptools.readthedocs.io/en/latest/setuptools.html#including-data-files
