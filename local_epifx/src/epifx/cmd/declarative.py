"""
Provides a declarative means of defining forecasts.

The purpose of this module is to allow users to define and run forecasting
simulations **without writing any Python code**, by instead defining all of
the necessary settings and parameters in a `TOML`_ file.
"""

import argparse
import logging

from .. import version


def common_parser(scenarios, config):
    """A command-line argument parser with common settings."""
    if scenarios:
        usage_str = '%(prog)s [options] scenario [scenario ...]'
    else:
        usage_str = '%(prog)s [options]'

    parser = argparse.ArgumentParser(usage=usage_str, add_help=False)

    h = parser.add_argument_group("Information")
    h.add_argument(
        '-h', '--help', action='help',
        help='Show this help message')
    h.add_argument("--version", action="version",
                   version="epifx {}".format(version.__version__),
                   help="Print the version information and exit.")

    og = parser.add_argument_group('Output options')
    log = og.add_mutually_exclusive_group()
    log.add_argument(
        '-d', '--debug', action='store_const', dest='loglevel',
        help='Enable debugging output',
        const=logging.DEBUG, default=logging.INFO)
    log.add_argument(
        '-q', '--quiet', action='store_const', dest='loglevel',
        help='Suppress logging output',
        const=logging.WARNING)

    if config:
        cg = parser.add_argument_group('Configuration settings')
        cg.add_argument(
            '-c', '--config', metavar='FILE', action='append',
            help=('The configuration file that defines the scenarios'))

    if scenarios:
        parser.add_argument(
            'scenario', nargs='*', metavar='scenario',
            help='Scenario identifier(s)')

    return parser
