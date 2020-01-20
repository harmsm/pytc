__description__ = \
"""
Classes describing the linear dependence of a free energy on osmolyte
concentration.
"""
__author__ = "Michael J. Harms"
__date__ = "2020-01-17"

from . import GlobalConnector

class LinearOsmolyte(GlobalConnector):
    """
    Returns the free energy for a reaction as a linear function of osmolyte.
    Requires experiment instances have .osm attribute.
    """

    param_guesses = {"m":0.0, "dG0":-5.0}
    required_data = ["osm"]

    def dG(self,experiment):
        """
        Return the change in free energy for a reaction as a function of
        osmolyte concentration.
        """

        return self.dG0 + self.m*experiment.osm

class LinearUrea(GlobalConnector):
    """
    Returns the free energy for a reaction as a linear function of urea.
    Requires experiment instances have .urea attribute.
    """

    param_guesses = {"m":0.0, "dG0":-5.0}
    required_data = ["urea"]

    def dG(self,experiment):
        """
        Return the change in free energy for a reaction as a function of
        urea concentration.
        """

        return self.dG0 + self.m*experiment.urea

class LinearGdm(GlobalConnector):
    """
    Returns the free energy for a reaction as a linear function of gdm.
    Requires experiment instances have .gdm attribute.
    """

    param_guesses = {"m":0.0, "dG0":-5.0}
    required_data = ["gdm"]

    def dG(self,experiment):
        """
        Return the change in free energy for a reaction as a function of
        gdm concentration.
        """

        return self.dG0 + self.m*experiment.gdm
