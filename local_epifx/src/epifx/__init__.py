"""Epidemic forecasts using disease surveillance data."""

import logging

from . import model
from . import det
from . import stoch
from . import obs
from . import summary
from . import select
from . import version

__package_name__ = u'epifx'
__author__ = u'Rob Moss'
__email__ = u'rgmoss@unimelb.edu.au'
__copyright__ = u'2014-2020, Rob Moss'
__license__ = u'BSD 3-Clause License'
__version__ = version.__version__


# Export classes from this module.
Model = model.Model

# Prevent an error message if the application does not configure logging.
log = logging.getLogger(__name__).addHandler(logging.NullHandler())
