__description__ = \
"""
A generic parameter that varies linearly.
"""
__author__ = "Michael J. Harms"
__date__ = "2020-01-17"

from . import GlobalConnector
import numpy as np

class GenericLinear(GlobalConnector):
    """
    Return a genric linear connector.
    """

    param_guesses = {"m":0.0, "b":0.0}
    required_data = []

    def __init__(self,name,x_name):
        """
        Initialize the class.

        name: name of the fitter.  will be pre-pended to all parameter names.
        x_name: name of the data with which the connector varies linearly
        """

        self._x_name = x_name
        self.required_data.append(x_name)

        super().__init__(name)

    def Y(self,experiment):
        """
        Return Y as a function of the values in x_name.
        """

        return self.b + self.m*experiment.__dict__[self._x_name]
