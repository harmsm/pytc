__description__ = \
"""

"""
__author__ = "Michael J. Harms"
__date__ = "2020-01-28"

from .base import GenericModel
from ..thermodynamics import two_state_melt

class TwoStateMelt(GenericModel):

    default_param_guesses = {"Ku":1e6,
                             "F_signal":1.0,
                             "U_signal":0.0}

    def __init__(self,df,Pt_column=None):

        # If no Pt column is specified, assume that Pt stays same across
        # melt
        if Pt_column is None:
            self._Pt = np.ones(len(df),dtype=np.float)

        super().__init__()

    @property
    def predicted(self):
        """
        Calculate the signal that would be observed.
        """

        Ku = self.param_values["Ku"]

        F, U = two_state_melt.species_conc(Ku,self._Pt)

        F_sig = F*self.param_values["F_signal"]
        U_sig = U*self.param_values["U_signal"]

        return F_sig + U_sig
