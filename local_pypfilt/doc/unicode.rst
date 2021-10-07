.. _unicode:

Unicode and byte strings
========================

The native ``str`` type is a Unicode string in
`Python 3 <https://docs.python.org/3/howto/unicode.html>`__.

.. tip::

    Software should only work with Unicode strings internally, decoding the
    input data as soon as possible and encoding the output only at the end
    (the "`Unicode sandwich <http://nedbatchelder.com/text/unipain.html>`__").

To that end, adhere to the following guidelines:

+ Use Unicode strings and Unicode literals everywhere.

+ If you have non-ASCII characters in a Python source file (e.g., in Unicode
  literals such as ``'α'``), you need to
  `declare the file encoding <https://www.python.org/dev/peps/pep-0263/>`__ at
  the top of the file:

  .. code-block:: python

     # -*- coding: utf-8 -*-

+ Encode Unicode text into UTF-8 when writing to disk:

  .. code-block:: python

     with open(filename, encoding='utf-8', mode='w') as f:
         f.write(unicode_string)

+ Decode UTF-8 bytes into Unicode text when reading from disk:

  .. code-block:: python

     with open(filename, encoding='utf-8') as f:
         unicode_lines = f.read().splitlines()

+ From NumPy 1.14 onward, functions such as
  `loadtxt <http://docs.scipy.org/doc/numpy-1.14.0/reference/generated/numpy.loadtxt.html>`__
  and
  `genfromtxt <http://docs.scipy.org/doc/numpy-1.14.0/reference/generated/numpy.genfromtxt.html>`__
  can handle files with arbitrary (Python-supported) text encoding:

  .. code-block:: python

     import numpy as np

     data = np.loadtxt(filename, encoding='utf-8', ...)

+ Use the ``'S'`` (byte string) data type when storing text in NumPy arrays.
  Encode Unicode text into UTF-8 when storing text, and decode UTF-8 bytes
  when reading text:

  .. code-block:: python

     >>> import numpy as np
     >>> xs = np.empty(3, dtype='S20')
     >>> xs[0] = 'abc'.encode('utf-8')
     >>> xs[1] = '« äëïöü »'.encode('utf-8')
     >>> xs[2] = 'ç'.encode('utf-8')
     >>> print(list(len(x) for x in xs))
     [3, 16, 2]
     >>> for x in xs:
     ...     print(x.decode('utf-8'))
     abc
     « äëïöü »
     ç

  .. note::

     There is also the option of using h5py's  `variable-length string type
     <http://docs.h5py.org/en/latest/special.html>`__ instead of ``'S'``.

+ NumPy has a Unicode data type (``'U'``), but it is
  `not supported by h5py <http://docs.h5py.org/en/latest/strings.html>`__ (and
  is platform-specific).

+ Note that h5py object names (i.e., groups and datasets) are
  `exclusively Unicode <http://docs.h5py.org/en/latest/strings.html>`__ and
  are stored as bytes, so byte strings will be used as-is and Unicode strings
  will be encoded using UTF-8.

+ Use Unicode strings and literals when encoding to and decoding from
  `JSON <http://blog.emacsos.com/unicode-in-python.html>`__:

  .. code-block:: python

     # Write UTF-8 bytes rather than '\uXXXX' escape sequences.
     with open(filename, encoding='utf-8', 'w') as f:
         json.dump(json_data, f, ensure_ascii=False)
