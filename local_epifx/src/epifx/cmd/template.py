"""Copy the example forecast settings files."""

import logging
import os.path

from . import settings
from .. import example


def parser():
    """Return the command-line argument parser for ``epifx-template``."""
    p = settings.common_parser(locns=False)

    tg = p.add_argument_group('Template options')
    tg.add_argument(
        'directory', nargs='?', default='.',
        help='The destination directory (default: "%(default)s")')

    return p


def main(args=None):
    """Copy the example forecast settings files."""
    p = parser()
    if args is None:
        args = vars(p.parse_args())
    else:
        args = vars(p.parse_args(args))
    logging.basicConfig(level=args['loglevel'])

    example.copy(args['directory'])
    www_dir = os.path.join(args['directory'], 'www')

    print("")
    print("NOTE: you must download 'd3' to view the JSON forecasts.")
    print("      https://github.com/d3/d3/releases/download/v3.5.17/d3.zip")
    print("")
    print("      Copy the file 'd3.min.js' to '{}'".format(www_dir))
    print("")
