
import numpy as np
import pandas as pd

def _species_conc(Kf,Ks,Kc,Pt,Ct):
    """
    Kf = D/U
    Ks = Ds/D
    Kc = DsC4/(Ds*C**4)
    Pt: total protein concentration
    Ct: total calcium concentration
    """

    # For construction of root polynomial
    p = np.zeros(6)
    p[0] = -Kc*Ks*Kf                     # C^5
    p[1] = (Ct*Kc*Ks*Kf - 4*Kc*Ks*Kf*Pt) # C^4
    p[2] = 0                             # C^3
    p[3] = 0                             # C^2
    p[4] = -(1 + Kf + Ks*Kf)
    p[5] = Ct*(1 + Kf + Ks*Kf)

    # Get root of C polynomial
    roots = np.roots(p)

    # Make sure we found a unique, real root that is between 0 and Ct
    mask = np.logical_and(np.logical_and(roots >= 0,roots <= Ct),np.isreal(roots))
    possible_roots = roots[mask]
    if len(possible_roots) == 0:
        #err = "no roots found\n"
        print("no roots found")
        return np.nan, np.nan, np.nan, np.nan, np.nan
        #raise RuntimeError(err)
    elif len(possible_roots) > 1:
        #err = "multiple roots found\n"
        print("multiple roots found")
        return np.nan, np.nan, np.nan, np.nan, np.nan
        #raise RuntimeError(err)
    else:
        pass

    # Determine C and U
    C = np.real(possible_roots[0])
    U = Pt/(1 + Kf + Ks*Kf + Kc*Ks*Kf*(C**4))
    D = Kf*U
    Ds = Ks*Kf*U
    DsC4 = Kc*Ks*Kf*U*(C**4)

    return U, D, Ds, DsC4, C

def species_conc(Kf,Ks,Kc,Pt,Ct):
    """
    Get species of U, D, Ds, and DsC4 given the equilibrium
    constants and the total protein and calcium concentrations.  This is the
    core, publically implementation of the model.  The equilibrium constants and
    concentrations may be floats or arrays, where all arrays have the same length.

    Kf: folding equilibrium constant Kf = D/U
    Ks: excited state equilibrium constant Ks = Ds/D
    Kc: calcium binding constant to all four sites Kc = DsC4/(Ds*C**4)
    Pt: total protein monomer concentration
    Ct: total calcium concentration

    returns U, D, Ds, DsC4
    """

    mismatch_error = False
    values = [Kf,Ks,Kc,Pt,Ct]
    is_array = [type(v) is pd.Series or type(v) is np.ndarray for v in values]

    # We have a series somewhere.  Make sure all series have the same length
    # and record this as series_length
    if sum(is_array) > 0:
        lengths = set([len(v) for i, v in enumerate(values) if is_array[i]])
        if len(lengths) != 1:
            mismatch_error = True
        else:
            series_length = list(lengths)[0]

    # No series.  series_length is 1.
    else:
        series_length = 1

    if mismatch_error:
        err = "Kf, Ks, Kc, Pt and Ct must either be float values or\n"
        err += "arrays.  Any arrays that are present must have the same length.\n"
        raise ValueError(err)

    # Go through values and either cast the series/array as an array or
    # repeat single values series_length times as arrays
    for i in range(len(values)):
        if is_array[i]:
            values[i] = np.array(values[i])
        else:
            values[i] = np.array([values[i] for _ in range(series_length)])

    # Pull values out of array
    Kf, Ks, Kc, Pt, Ct = values

    U, D, Ds, DsC4, C = [], [], [], [], []
    for i in range(len(Pt)):

        concs = _species_conc(Kf[i],Ks[i],Kc[i],Pt[i],Ct[i])

        U.append(concs[0])
        D.append(concs[1])
        Ds.append(concs[2])
        DsC4.append(concs[3])
        C.append(concs[4])

    U = np.array(U)
    D = np.array(D)
    Ds = np.array(Ds)
    DsC4 = np.array(DsC4)
    C = np.array(C)

    return U, D, Ds, DsC4, C
