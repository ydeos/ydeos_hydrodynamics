# coding: utf-8

r"""Domain validity checks."""

import logging
import warnings
from ydeos_hydrodynamics.memoize import memoize


logger = logging.getLogger(__name__)


@memoize
def check_domain(domain: str,
                 quantity: str,
                 lower: float,
                 upper: float,
                 value: float) -> bool:
    r"""Check a value is within the lower to upper range.

    IF it isn't, report using the domain and the quantity

    """
    if lower <= value <= upper:
        return True
    msg = f"{quantity} [{value}] is outside the %s range [{lower} - {upper}]"
    warnings.warn(msg)
    logger.warning(msg)
    return False
