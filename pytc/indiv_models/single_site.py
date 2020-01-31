__description__ = \
"""
base.py

Base class for all model description.
"""
__author__ = "Michael J. Harms"
__date__ = "2016-06-22"

from .base import PytcModel
from ..thermodynamics import single_site_binding

class SingleSite(PytcModel):

    default_param_guesses = {"Kd":1e-6,
                             "P_signal":0.0,
                             "X_signal":0.0,
                             "PX_signal":1.0}

    def __init__(self,df,Pt_column="Pt",Xt_column="Xt"):

        self._Pt = df[Pt_column]
        self._Xt = df[Xt_column]

        super().__init__()

    @property
    def predicted(self):
        """
        Calculate the signal that would be observed.
        """

        Kd = self.param_values["Kd"]

        P, X, PX = single_site_binding.species_conc(Kd,self._Pt,self._Xt)

        P_sig = P*self.param_values["P_signal"]
        X_sig = X*self.param_values["X_signal"]
        PX_sig = PX*self.param_values["PX_signal"]

        return P_sig + X_sig + PX_sig
