"""Compare HDF5 simulation files for identical outputs and metadata."""

import h5py
import logging
import numpy as np

from . import settings


def arrays(f1, f2, path, verbose=True):
    """Compare two HDF5 arrays for equality.

    :param f1: The first HDF5 file.
    :param f2: The second HDF5 file.
    :param path: The path of the array dataset.
    :param verbose: Whether to report successful matches.

    :return: ``True`` if the arrays are equal, otherwise ``False``.
    """
    try:
        arr1 = f1[path][()]
    except KeyError:
        arr1 = None

    try:
        arr2 = f2[path][()]
    except KeyError:
        arr2 = None

    if arr1 is None and arr2 is None:
        if verbose:
            print("Both missing {}".format(path))
            return True
    elif arr1 is None:
        print("Missing {} from {}".format(path, f1.filename))
        return False
    elif arr2 is None:
        print("Missing {} from {}".format(path, f2.filename))
        return False

    if arr1.shape != arr2.shape:
        print("Different shape for {}: {} vs {}".format(path, arr1.shape,
                                                        arr2.shape))
        return False
    elif not np.array_equal(arr1, arr2):
        print("Different values for {}".format(path))
        return False
    elif verbose:
        print("Match for {}".format(path))

    return True


def groups(f1, f2, path, verbose=True):
    """Compare two HDF5 groups for equality.

    :param f1: The first HDF5 file.
    :param f2: The second HDF5 file.
    :param path: The path of the group.
    :param verbose: Whether to report successful matches.

    :return: ``True`` if the groups are equal, otherwise ``False``.
    """
    try:
        grp1 = f1[path]
    except KeyError:
        grp1 = None

    try:
        grp2 = f2[path]
    except KeyError:
        grp2 = None

    if grp1 is None and grp2 is None:
        if verbose:
            print("Both missing {}".format(path))
            return True
    elif grp1 is None:
        print("Missing {} from {}".format(path, f1.filename))
        return False
    elif grp2 is None:
        print("Missing {} from {}".format(path, f2.filename))
        return False

    keys1, keys2 = set(grp1.keys()), set(grp2.keys())
    diff1, diff2 = keys1.difference(keys2), keys2.difference(keys1)

    same = True

    if diff1:
        print("Missing keys for {} in {}:\n    ".format(path, f2.filename),
              [x for x in diff1])
        same = False
    if diff2:
        print("Missing keys for {} in {}:\n    ".format(path, f1.filename),
              [x for x in diff2])
        same = False

    for k in keys1.intersection(keys2):
        val1, val2 = grp1[k], grp2[k]
        k_path = "{}/{}".format(path, k)
        # Ignore the command line, since it can differ without affecting the
        # simulation (e.g., by changing the order of any optional arguments).
        if k_path == '/meta/sim/cmdline':
            continue
        if isinstance(val1, h5py.Dataset) and isinstance(val2, h5py.Dataset):
            if val1.shape == () and val2.shape == ():
                if val1.dtype != val2.dtype:
                    print("Different types for {}: '{}' vs '{}'".format(
                        k_path, val1.dtype, val2.dtype))
                    same = False

                val1, val2 = val1[()], val2[()]
                if val1 != val2:
                    print("Different values for {}: '{}' vs '{}'".format(
                        k_path, val1, val2))
                    same = False
                elif verbose:
                    print("Match for {}".format(k_path))
            else:
                same = arrays(f1, f2, k_path, verbose=verbose) and same
        elif isinstance(val1, h5py.Group) and isinstance(val2, h5py.Group):
            same = groups(f1, f2, k_path, verbose=verbose) and same
        else:
            print("Different types for {}: {} vs {}".format(path, type(val1),
                                                            type(val2)))
            same = False

    return same


def files(path1, path2, verbose=True, examine=None):
    """Compare two HDF5 files for identical simulation outputs.

    :param path1: The filename of the first HDF5 file.
    :param path2: The filename of the second HDF5 file.
    :param verbose: Whether to report successful matches.
    :param examine: The data groups to examine for equality; the default is to
        examine the simulation outputs (``'/data'``) and ignore the associated
        metadata (``'/meta'``).

    :return: ``True`` if the files contain identical simulation outputs,
        otherwise ``False``.
    """
    same = True
    if examine is None:
        # Default to comparing only the simulation ouputs.
        examine = ['/data']
    with h5py.File(path1, 'r') as f1, h5py.File(path2, 'r') as f2:
        for key in examine:
            # Compare the specified data group.
            same = groups(f1, f2, key, verbose=verbose) and same

    return same


def parser():
    """The command-line argument parser for ``epifx.cmp.main``."""
    p = settings.common_parser()
    p.description = "Compare HDF5 files for identical simulation outputs"

    # Find the previously-created "Output options" group, if it exists.
    gs = [g for g in p._action_groups if g.title == "Output options"]
    if len(gs) == 1:
        g = gs[0]
    else:
        g = p.add_argument_group("Output options")

    g.add_argument("-v", "--verbose", action="store_true",
                   help="Report identical groups and datasets.")
    g.add_argument("-m", "--meta-data", action="store_const",
                   const=["/data", "/meta"], dest="examine",
                   help="Compare output data *and* metadata.")

    f = p.add_argument_group("Input files")
    f.add_argument('file1', metavar='file1.hdf5',
                   help='The simulation files to compare.')
    f.add_argument('file2', metavar='file2.hdf5',
                   help='The simulation files to compare.')

    return p


def main(args=None):
    p = parser()
    if args is None:
        args = vars(p.parse_args())
    else:
        args = vars(p.parse_args(args))
    logging.basicConfig(level=args['loglevel'])

    same = files(args['file1'], args['file2'], verbose=args['verbose'],
                 examine=args['examine'])

    if same:
        print("Files match")
        return 0
    else:
        print("Files do not match")
        return 10
