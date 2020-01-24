__description__ = \
"""
base.py

Base class for other itc model description.
"""
__author__ = "Michael J. Harms"
__date__ = "2016-06-22"

import numpy as np
from ..base import PytcModel

class ITCModel(PytcModel):
    """
    Base class from which all ITC models should be sub-classed.
    """

    def __init__(self,
                 S_cell=100e-6,S_syringe=0.0,
                 T_cell=0.0,   T_syringe=1000e-6,
                 cell_volume=300.0,
                 shot_volumes=[2.5 for i in range(30)]):

        """
        S_cell: stationary concentration in cell in M
        S_syringe: stationary concentration in syringe in M
        T_cell: titrant concentration cell in M
        T_syringe: titrant concentration syringe in M
        cell_volume: cell volume, in uL
        shot_volumes: list of shot volumes, in uL.
        """

        self._S_cell = S_cell
        self._S_syringe = S_syringe

        self._T_cell = T_cell
        self._T_syringe = T_syringe

        self._cell_volume = cell_volume
        self._shot_volumes = np.array(shot_volumes)

        # Determine the concentration of all of the species across the titration
        self._S_conc = self._titrate_species(self._S_cell,self._S_syringe)
        self._T_conc = self._titrate_species(self._T_cell,self._T_syringe)

        # Add dilution parameters
        self._initialize_param(param_names=["dilution_heat","dilution_intercept"],
                               param_guesses=[0.0,0.0])

    @property
    def predicted(self):
        return np.array(())

    # --------------------------------------------------------------------------

    def _titrate_species(self,cell_conc,syringe_conc):
        """
        Determine the concentrations of stationary and titrant species in the
        cell given a set of titration shots and initial concentrations of both
        the stationary and titrant species.

        Does two independent calculations and adds them.  First, it calculates
        the concentration change due to injection (injected).  It then treats
        the dilution of the stuff initially in hte cell (diluted).  The sum of
        these two groups is the total concentration of whatever was titrated.
        The shot_ratio_product method is described on p. 134 of Freire et al.
        (2009) Meth Enzymology 455:127-155
        """

        volume = np.zeros(len(self._shot_volumes)+1)
        out_conc = np.zeros(len(self._shot_volumes)+1)

        volume[0] = self._cell_volume
        out_conc[0] = cell_conc

        shot_ratio = (1 - self._shot_volumes/self._cell_volume)
        for i in range(len(self._shot_volumes)):

            shot_ratio_prod = np.prod(shot_ratio[:(i+1)])
            injected = syringe_conc*(1 - shot_ratio_prod)
            diluted = cell_conc*shot_ratio_prod

            out_conc[i+1] = injected + diluted

        return out_conc

    @property
    def mole_ratio(self):
        """
        Molar ratio of titrant to stationary species.  If not yet initialized,
        send return empty array.
        """

        try:
            return self._T_conc[1:]/self._S_conc[1:]
        except AttributeError:
            return np.array([],dtype=float)

    @property
    def dilution_heats(self):
        """
        Return the heat of dilution.
        """

        return self._T_conc[1:]*self.param_values["dilution_heat"] + self.param_values["dilution_intercept"]
