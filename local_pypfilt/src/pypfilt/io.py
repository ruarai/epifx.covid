"""Read data from external sources."""

import datetime
import numpy as np


def read_table(path, columns, comment='#', encoding='utf-8'):
    """
    Read data from a space-delimited text file with column headers defined in
    the first non-comment line.

    :param path: The path to the data file.
    :param columns: The columns to read from the data file, represented as a
        sequence of ``(name, type)`` tuples where ``type`` must be a NumPy
        `scalar type <https://numpy.org/doc/stable/reference/arrays.scalars.html>`__,
        or ``(name, type, converter)`` tuples where ``converter`` is a
        function that converts the column string into the desired value.
    :param comment: The characters, or list of characters, that indicate the
        start of a single-line comment.
    :param encoding: The name of the encoding used to decode the file content.

    :Examples:

    >>> from pypfilt.io import date_column, read_table
    >>> import numpy as np
    >>> import datetime
    >>> path = "input_data.ssv"
    >>> with open(path, 'w') as f:
    ...    _ = f.write('date value\\n')
    ...    _ = f.write('2020-01-01 1\\n')
    ...    _ = f.write('2020-01-02 3\\n')
    ...    _ = f.write('2020-01-03 5\\n')
    >>> columns = [date_column('date'), ('value', np.int_)]
    >>> data = read_table(path, columns)
    >>> isinstance(data['date'][0], datetime.datetime)
    True
    >>> observations = [{'date': row['date'], 'value': row['value']}
    ...                 for row in data]
    """

    # Read in the column names, and find the line where the data begins.
    skip_lines = 1
    with open(path, encoding=encoding) as f:
        file_cols = f.readline().strip().split()
        while len(file_cols) == 0 or file_cols[0].startswith(comment):
            file_cols = f.readline().strip().split()
            skip_lines += 1

    # Ensure all of the required columns are defined.
    need_columns = [col_tuple[0] for col_tuple in columns]
    for column_name in need_columns:
        if column_name not in file_cols:
            raise ValueError('Column "{}" not found in {}'.format(
                column_name, path))

    # Construct the list of column types and associate column conversion
    # functions with the index of that column in the file.
    converters = {}
    column_dtypes = []
    for ix, col_tuple in enumerate(columns):
        if len(col_tuple) == 2:
            column_dtypes.append(col_tuple)
        elif len(col_tuple) == 3:
            column_dtypes.append(col_tuple[:2])
            column_ix = file_cols.index(col_tuple[0])
            converters[column_ix] = col_tuple[2]

    # Determine the index of each required column in the file.
    read_columns = [file_cols.index(name) for name in need_columns]

    tbl = np.loadtxt(path, encoding=encoding,
                     skiprows=skip_lines, dtype=column_dtypes,
                     converters=converters, usecols=read_columns)

    # If the table only contains a single row, it will be represented as a
    # scalar value. So we need to convert it to an array with at least one
    # dimension.
    return np.atleast_1d(tbl)


def date_column(name, fmt='%Y-%m-%d'):
    """
    Return a ``(name, type, converter)`` tuple that can be used with
    :func:`read_table` to convert a column into ``datetime.datetime`` values.

    .. note::

       Where dates are used for observation times, they should be represented
       as ``datetime.datetime`` values, not as ``datetime.date`` values.
       This is why this function returns a converter that returns
       ``datetime.datetime`` values.

    :param str name: The column name in the data file.
    :param str fmt: The date format used to parse the column values.
    """
    return (name, np.object_,
            lambda s: datetime.datetime.strptime(s, fmt))


def datetime_column(name, fmt='%Y-%m-%dT%H:%M:%S'):
    """
    Return a ``(name, type, converter)`` tuple that can be used with
    :func:`read_table` to convert a column into ``datetime.datetime`` values.

    :param str name: The column name in the data file.
    :param str fmt: The datetime format used to parse the column values.
    """
    return (name, np.object_,
            lambda s: datetime.datetime.strptime(s, fmt))


def read_lookup_table(path, time, dtype='f8', comment='#', encoding='utf-8'):
    """
    Read time-indexed data from a space-delimited text file with column
    headers defined in the first non-comment line.

    :param path: The path to the data file.
    :param pypfilt.time.Time time: The time scale.
    :param dtype: The type of the lookup values.
    :param comment: The characters, or list of characters, that indicate the
        start of a single-line comment.
    :param encoding: The name of the encoding used to decode the file content.

    :Examples:

    >>> from pypfilt.io import read_lookup_table, lookup
    >>> from pypfilt.time import Datetime
    >>> import datetime
    >>> path = "input_data.ssv"
    >>> with open(path, 'w') as f:
    ...    _ = f.write('date value1 value2 value3\\n')
    ...    _ = f.write('2020-01-01 1.0 1.5 2.0\\n')
    >>> time = Datetime()
    >>> table = read_lookup_table(path, time)
    >>> isinstance(table['date'][0], datetime.datetime)
    True
    >>> when = datetime.datetime(2020, 1, 1)
    >>> values = lookup(table, when)
    >>> len(values.shape) == 1
    True
    >>> all(isinstance(value, float) for value in values)
    True
    """

    # Read in the column names, and find the line where the data begins.
    skip_lines = 1
    with open(path, encoding='utf-8') as f:
        cols = f.readline().strip().split()
        while len(cols) == 0 or cols[0].startswith(comment):
            cols = f.readline().strip().split()
            skip_lines += 1

    columns = [(name, dtype) for name in cols[1:]]
    columns.insert(0, time.column(cols[0]))

    tbl = read_table(path, columns)

    # NOTE: rename the first column to 'date', so that the cache can identify
    # this as a time-indexed table.
    col_names = list(tbl.dtype.names)
    col_names[0] = 'date'
    tbl.dtype.names = col_names

    # If the table only contains a single row, it will be represented as a
    # scalar value. So we need to convert it to an array with at least one
    # dimension.
    tbl = np.atleast_1d(tbl)

    if len(tbl) == 0:
        raise ValueError("File '{}' contains no rows".format(path))

    # Count all columns except for 'date'.
    num_value_cols = len(tbl[0]) - 1

    # Transform this table with rows (date, value1, value2, ...) into a table
    # with rows (date, [values]).
    new_dtype = [time.column('date')[:2],
                 ('value', dtype, (num_value_cols,))]
    row_pairs = [(row[0], tuple(row)[1:]) for row in tbl]
    value_tbl = np.array(row_pairs, dtype=new_dtype)

    return value_tbl


def lookup_values_count(lookup_table):
    """
    Return the number of value columns in a lookup table.
    """
    # NOTE: ignore the first column, which contains the lookup time.
    return lookup_table.dtype['value'].shape[0]


def lookup(lookup_table, when):
    """
    Return the values associated with a specific time.
    """
    time_col = lookup_table.dtype.names[0]
    ixs = np.where(lookup_table[time_col] <= when)[0]
    if len(ixs) == 0:
        # No match, default to the earliest value.
        most_recent_row = 0
    else:
        # Otherwise use the most recent match.
        most_recent_row = ixs[-1]
    return lookup_table[most_recent_row]['value']


class Lookup:
    def __init__(self, lookup_table):
        self.__table = lookup_table

    def value_count(self):
        return lookup_values_count(self.__table)

    def lookup(self, when):
        return lookup(self.__table, when)
