__description__ = \
"""
experiments.py

Classes for loading experimental data and associating those data with a
model.

Units:
    Volumes are in microliters
    Temperatures are in Kelvin
    Concentrations are in molar
    Energy is `units`, where `units` is specified when instantiating the
    PytcExperiment class.  It must be a in the AVAIL_UNITS dictionary.
"""
__author__ = "Michael J. Harms"
__date__ = "2016-06-22"

import random, string, os
import numpy as np

class PytcExperiment:
    """
    Class that holds an experimental measurement and a model that describes it.
    """

    AVAIL_UNITS = {"cal/mol":1.9872036,
                   "kcal/mol":0.0019872036,
                   "J/mol":8.3144598,
                   "kJ/mol":0.0083144598}

    def __init__(self,data_file,model,units="cal/mol",
                 uncertainty=0.1,**model_kwargs):
        """

        Parameters
        ----------

        data_file: string
            file containing experimental results
        model: PytcModel subclass instance
            PytcModel subclass to use for modeling
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
                        stored in the dh_file.
        """

        self.dh_file = dh_file
        self._shot_start = shot_start

        # Deal with units
        self._units = units
        try:
            self._R = self.AVAIL_UNITS[self._units]
        except KeyError:
            err = "units must be one of:\n"
            for k in self.AVAIL_UNITS.keys():
                err += "    {}\n".format(k)
            err += "\n"

            raise ValueError(err)

        # For numerical reasons, there should always be *some* uncertainty
        self._uncertainty = uncertainty
        if self._uncertainty == 0.0:
            self._uncertainty = 1e-12

        # Load in heats
        extension = self.dh_file.split(".")[-1]
        self._read_heats_file()

        # Initialize model using information read from heats file
        self._model = model(S_cell=self.stationary_cell_conc,
                            T_syringe=self.titrant_syringe_conc,
                            cell_volume=self.cell_volume,
                            shot_volumes=self._shots,**model_kwargs)

        r = "".join([random.choice(string.ascii_letters) for i in range(20)])
        self._experiment_id = "{}_{}".format(self.dh_file,r)


    def _read_heats_file(self):
        """
        Dummy heat reading file.
        """

        pass

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
    def param_values(self):
        """
        Values of fit parameters.
        """

        return self._model.param_values

    @property
    def param_stdevs(self):
        """
        Standard deviations on fit parameters.
        """

        return self._model.param_stdevs

    @property
    def param_ninetyfives(self):
        """
        95% confidence intervals on fit parmeters.
        """

        return self._model.param_ninetyfives

    @property
    def model(self):
        """
        Fitting model.
        """

        return self._model

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
    def obs(self):
        """
        Return experimental observable.
        """
        return self._obs[:]

    @heats.setter
    def obs(self,obs):
        """
        Set the obs.
        """

        self._obs[:] = obs[:]

    @property
    def obs_stdev(self):
        """
        Standard deviation on the uncertainty of the observable.
        """

        return self._obs_stdev[:]

    @obs_stdev.setter
    def obs_stdev(self,obs_stdev):
        """
        Set the standard deviation on the uncertainty of the observable.
        """

        self._obs_stdev[:] = obs_stdev[:]

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

    @property
    def experiment_id(self):
        """
        Return a unique experimental id.
        """

        return self._experiment_id

    @property
    def units(self):
        """
        Units for file.
        """

        return self._units

    @units.setter
    def units(self,units):
        """
        Change the units.
        """

        # Deal with units
        self._units = units
        try:
            self._R = self.AVAIL_UNITS[self._units]
        except KeyError:
            err = "units must be one of:\n"
            for k in self.AVAIL_UNITS.keys():
                err += "    {}\n".format(k)
            err += "\n"

            raise ValueError(err)

    @property
    def R(self):
        """
        Experiment gas constant.
        """

        return self._R
