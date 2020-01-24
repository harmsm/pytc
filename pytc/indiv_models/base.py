__description__ = \
"""
base.py

Base class for all model description.
"""
__author__ = "Michael J. Harms"
__date__ = "2016-06-22"

import numpy as np
from .. import fit_param

class PytcModel:
    """
    Base class from which all pytc models should be sub-classed.
    """

    def __init__(self):
        """
        """

        # Add dilution parameters
        self._initialize_param()

    def _initialize_param(self,param_names=None,param_guesses=None):
        """
        Initialize the parameters.
        """

        self._params = {}

        if param_names == None:
            param_names = []
        if param_guesses == None:
            param_guesses = []

        if len(param_names) != len(param_guesses):
            err = "parameter names and parameter guesses must have the same\n"
            err += "length.\n"
            raise ValueError(err)

        # Grab parameter names and guesses from self.default_param_guesses (if
        # specified).  Anything specified in _initialize_param arguments will
        # be taken rather than these values.
        try:

            for p in self.default_param_guesses.keys():
                if p not in param_names:
                    param_names.append(p)
                    param_guesses.append(self.default_param_guesses[p])
        except AttributeError:
            pass

        for i, p in enumerate(param_names):
            self._params[p] = fit_param.FitParameter(p,guess=param_guesses[i])

        self._param_names = param_names[:]
        self._param_names.sort()


    @property
    def predicted(self):
        """
        Predicted observations given model. 
        """
        return self._model.predicted

    # -------------------------------------------------------------------------
    # parameter names

    @property
    def param_names(self):
        """
        The parameters for the model.
        """

        return self._param_names

    # -------------------------------------------------------------------------
    # parameter objects

    @property
    def parameters(self):
        """
        Return FitParam objects associated with the model.
        """

        return self._params

    # -------------------------------------------------------------------------
    # parameter values

    @property
    def param_values(self):
        """
        Values for each parameter in the model.
        """

        return dict([(p,self._params[p].value) for p in self._param_names])


    def update_values(self,param_values):
        """
        Update parameter values for fit. param_values is a dictionary with
        with some number of parameter names.
        """

        for p in param_values.keys():
            self._params[p].value = param_values[p]

    # -------------------------------------------------------------------------
    # parameter stdev

    @property
    def param_stdevs(self):
        """
        Standard deviation for each parameter in the model.
        """

        return dict([(p,self._params[p].stdev) for p in self._param_names])


    def update_stdevs(self,param_stdevs):
        """
        Update parameter stdev for fit. param_stdevs is a dictionary with
        with some number of parameter names.
        """

        for p in param_stdevs.keys():
            self._params[p].stdev = param_stdevs[p]

    # -------------------------------------------------------------------------
    # parameter ninetyfive

    @property
    def param_ninetyfives(self):
        """
        95% confidence intervals for each parameter in the model.
        """

        return dict([(p,self._params[p].ninetyfive) for p in self._param_names])


    def update_ninetyfives(self,param_ninetyfives):
        """
        Update parameter 95% for fit. param_ninetyfives is a dictionary with
        with some number of parameter names.
        """

        for p in param_ninetyfives.keys():
            self._params[p].ninetyfive = param_ninetyfives[p]

    # -------------------------------------------------------------------------
    # parameter guesses

    @property
    def param_guesses(self):
        """
        Guesses for each parameter in the model.
        """

        return dict([(p,self._params[p].guess) for p in self._param_names])

    def update_guesses(self,param_guesses):
        """
        Update parameter guesses for fit. param_guesses is a dictionary with
        with some number of parameter names.
        """

        for p in param_guesses.keys():
            self._params[p].guess = param_guesses[p]

    # -------------------------------------------------------------------------
    # parameter ranges

    @property
    def param_guess_ranges(self):
        """
        Return parameter ranges.
        """

        return dict([(p,self._params[p].guess_range) for p in self._param_names])

    def update_guess_ranges(self,param_ranges):
        """
        Update parameter ranges.  param_ranges is a dictionary of paramters
        keyed to two-entry lists/tuples or ranges.
        """

        for p in param_ranges.keys():
            self._params[p].guess_range = param_ranges[p]


    # -------------------------------------------------------------------------
    # fixed parameters

    @property
    def fixed_param(self):
        """
        Return the fixed parameters.
        """

        return dict([(p,self._params[p].fixed) for p in self._param_names])

    def update_fixed(self,fixed_param):
        """
        Fix parameters.  fixed_param is a dictionary of parameters keyed to their
        fixed values.  If the value is None, the parameter is removed from the
        fixed parameters dictionary and will float.
        """

        for p in fixed_param.keys():

            if fixed_param[p] == None:
                self._params[p].fixed = False
            else:
                self._params[p].fixed = True
                self._params[p].value = fixed_param[p]


    # -------------------------------------------------------------------------
    # parameter bounds

    @property
    def bounds(self):
        """
        Return parameter bounds.
        """

        return dict([(p,self._params[p].bounds) for p in self._param_names])

    def update_bounds(self,bounds):
        """
        Update parameter bounds.  bounds is a dictionary of paramters
        keyed to two-entry lists/tuples or ranges.
        """

        for p in bounds.keys():
            self._params[p].bounds = bounds[p]

    # -------------------------------------------------------------------------
    # parameter aliases

    @property
    def param_aliases(self):
        """
        Return parameter aliases.
        """

        return dict([(p,self._params[p].alias) for p in self._param_names
                     if self._params[p].alias != None])

    def update_aliases(self,param_alias):
        """
        Update parameter aliases.  param_alias is a dictionary of parameters keyed
        to their aliases (used by the global fit).  If the value is None, the parameter
        alias is removed.
        """

        for p in param_alias.keys():
            self._params[p].alias = param_alias[p]
