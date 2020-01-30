__description__ = \
"""
Generic log transformation of a parameter.
"""
__author__ = "Michael J. Harms"
__date__ = "2020-01-17"

from . import GlobalConnector
import numpy as np

class LogTransform(GlobalConnector):
    """
    Return a genric linear connector.
    """

    param_guesses = {"logK":0.0}

    def exp_logK(self,experiment):
        """
        Return exp(logK)
        """

        return np.exp(self.logK)

    def ten_logK(self,experiment):
        """
        Return 10**logK.
        """

        return 10**(self.logK)
