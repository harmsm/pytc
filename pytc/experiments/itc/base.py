__description__ = \
"""
experiments.py

Classes for loading experimental ITC data and associating those data with a
model.

Units:
    Volumes are in microliters
    Temperatures are in Kelvin
    Concentrations are in molar
    Energy is `units`, where `units` is specified when instantiating the
    ITCExperiment class.  It must be a in the AVAIL_UNITS dictionary.
"""
__author__ = "Michael J. Harms"
__date__ = "2016-06-22"

import random, string, os
import numpy as np

from ..base import PytcExperiment

class BaseITCExperiment(PytcExperiment):
    """
    Class that holds an experimental ITC measurement and a model that describes it.
    """

    AVAIL_UNITS = {"cal/mol":1.9872036,
                   "kcal/mol":0.0019872036,
                   "J/mol":8.3144598,
                   "kJ/mol":0.0083144598}

    def __init__(self,data_file,model,shot_start=1,units="cal/mol",
                 uncertainty=0.1,**model_kwargs):
        """

        Parameters
        ----------

        data_file: string
            integrated heats file written out by origin software.
        model: ITCModel subclass instance
            ITCModel subclass to use for modeling
        shot_start: int
            what shot to use as the first real point.  Shots start at 0, so
            default=1 discards first point.
        units : string
            file units ("cal/mol","kcal/mol","J/mol","kJ/mol")
        uncertainty : float > 0.0
            uncertainty in integrated heats (set to same for all shots, unless
            specified in something like NITPIC output file).

        **model_kwargs: any keyword arguments to pass to the model.  Any
                        keywords passed here will override whatever is
                        stored in the data_file.
        """


        # ITC specific information
        self._shot_start = shot_start

        # Run parent init
        super().__init__(data_file,model,units,uncertainty,**model_kwargs)


    def _initialize_model(self,**model_kwargs):

        # Initialize model using information read from heats file
        self._model = self._model_function(S_cell=self.stationary_cell_conc,
                                           T_syringe=self.titrant_syringe_conc,
                                           cell_volume=self.cell_volume,
                                           shot_volumes=self._shots,
                                           **model_kwargs)

    @property
    def dQ(self):
        """
        Return heats calculated by the model with parameters defined in params
        dictionary.
        """

        if len(self._model.dQ) == 0:
            return np.array(())

        return self._model.dQ[self._shot_start:]

    @property
    def dilution_heats(self):
        """
        Return dilution heats calculated by the model with parameters defined
        in params dictionary.
        """

        if len(self._model.dilution_heats) == 0:
            return np.array(())

        return self._model.dilution_heats[self._shot_start:]


    @property
    def shot_start(self):
        """
        Starting shot to use.
        """

        return self._shot_start

    @shot_start.setter
    def shot_start(self,value):
        """
        Change starting shot.
        """

        self._shot_start = value

    @property
    def heats(self):
        """
        Return experimental heats.
        """
        return self._heats[self._shot_start:]

    @heats.setter
    def heats(self,heats):
        """
        Set the heats.
        """

        self._heats[self._shot_start:] = heats[:]

    @property
    def heats_stdev(self):
        """
        Standard deviation on the uncertainty of the heat.
        """

        return self._heats_stdev[self._shot_start:]

    @heats_stdev.setter
    def heats_stdev(self,heats_stdev):
        """
        Set the standard deviation on the uncertainty of the heat.
        """

        self._heats_stdev[self._shot_start:] = heats_stdev[:]

    @property
    def mol_injected(self):
        """
        Return the mols injected over shots.
        """

        # uL * mol/L * L/1e6 uL -> mol
        return self._shots[self._shot_start:]*self.titrant_syringe_conc*1e-6

    @property
    def mole_ratio(self):
        """
        Return the mole ratio of titrant to stationary.
        """
        return self._model.mole_ratio[self._shot_start:]
