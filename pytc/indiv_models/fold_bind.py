__description__ = \
"""
Linked folding, dimerization, and binding model.  Written to describe behavior
of the protein human S100A9.
"""
__author__ = "Michael J. Harms"
__date__ = "2020-01-28"

from .base import PytcModel
from ..thermodynamics import fold_bind

import numpy as np

class FoldBind(PytcModel):

    default_param_guesses = {"Kf":1e6,
                             "Ks":1e6,
                             "Kc":1e12,
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
        Ks = self.param_values["Ks"]
        Kc = self.param_values["Kc"]

        Pt = self._Pt
        Ct = self._Ct

        U, D, Ds, DsC4, C = fold_bind.species_conc(Kf,Ks,Kc,Pt,Ct)

        Ds_sig = (Ds + DsC4)*self.param_values["Ds_signal"]
        D_sig = D*self.param_values["D_signal"]
        U_sig = U*self.param_values["U_signal"]

        return (Ds_sig + D_sig + U_sig)/Pt