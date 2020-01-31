__description__ = \
"""
Linked folding, dimerization, and binding model.  Written to describe behavior
of the protein human S100A9.
"""
__author__ = "Michael J. Harms"
__date__ = "2020-01-28"

from .base import PytcModel
from ..thermodynamics import fold_bind_dimerize

import numpy as np

class FoldBindDimerize(PytcModel):

    default_param_guesses = {"Kf":1e6,
                             "Kd":1e6,
                             "Ks":0.1,
                             "K1":1e6,
                             "K2":1e6,
                             "Ds_signal":1.0,
                             "D_signal":1.0,
                             "U_signal":0.0}

    def __init__(self,df,Pt_column="Pt",Ct_column="Ct"):

        self._Pt = df[Pt_column]
        self._Ct = df[Ct_column]

        super().__init__()

    @property
    def predicted(self):
        """
        Calculate the signal that would be observed.  This is assumes a
        molar signal (mean-molar ellipticity or the like), not an absolute
        signal. This matters if you have different experiments at different
        protein concentrations.
        """

        Kf = self.param_values["Kf"]
        Kd = self.param_values["Kd"]
        Ks = self.param_values["Ks"]
        K1 = self.param_values["K1"]
        K2 = self.param_values["K2"]


        Pt = self._Pt
        Ct = self._Ct

        U, M, D, Ds, DsC1, DsC2, DsC3, DsC4, C = fold_bind_dimerize.species_conc(Kf,Kd,Ks,K1,K2,Pt,Ct)

        Ds_sig = (Ds + DsC1 + DsC2 + DsC3 + DsC4)*self.param_values["Ds_signal"]
        D_sig = D*self.param_values["D_signal"]
        U_sig = U*self.param_values["U_signal"]

        return (Ds_sig + D_sig + U_sig)/Pt
