__description__ = \
"""
Global models for fitting collections of experiments.
"""
__author__ = "Michael J. Harms"
__date__ = "2016-12-19"
__all__ = []

from .base import GlobalConnector

from .num_protons import NumProtons
from .vant_hoff import VantHoff
from .vant_hoff_extended import VantHoffExtended
from .linear_osmolyte import LinearOsmolyte, LinearUrea, LinearGdm
from .generic_linear import GenericLinear
