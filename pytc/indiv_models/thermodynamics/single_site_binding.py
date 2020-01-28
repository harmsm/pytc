
import numpy as np

def species_conc(Kd,Pt,Xt):
    """
    Single-site binding.

    Kd: equilibrium constant for dissociation
    Pt: total protein concentration
    Xt: total ligand concentration

    Returns:
    P: free protein concentration
    X: free ligand concentration
    PX: ligand bound concentration
    """

    b = Pt + Xt + Kd

    PX = (b - np.sqrt((b)**2 - 4*Pt*Xt))/2
    P = Pt - PX
    X = Xt - PX

    return P, X, PX
