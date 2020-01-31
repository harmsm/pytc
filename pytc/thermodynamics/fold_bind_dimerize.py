__description__ = \
"""
Linked folding, dimerization, and binding model.  Written to describe behavior
of the protein human S100A9.
"""
__author__ = "Michael J. Harms (harmsm@gmail.com)"
__date__ = "2020-01-15"

import numpy as np
import pandas as pd
from scipy import optimize

def _is_tiny(x,Pt,tol=1e-7):
    """
    Determine whether x makes a measurable contribution to Pt.
    """

    if x < Pt*tol:
        return True
    else:
        return False


def _calc_prot_conc(Kf,Kd,Ks,K1,K2,Pt,C):
    """
    Calculate the concentration of all protein species given the equilibrium
    constants, total protein concentration (monomer), and free calcium
    concentration.
    """

    # If we have no protein, we have no subspecies
    if Pt == 0:
        return 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0

    # Binding polynomial for dimeric species
    Dbar = 1 + Ks*(1 + 2*K1*C + K1*K1*C*C + 2*K1*K1*K2*C*C*C + K1*K1*K2*K2*C*C*C*C)

    # If we really have no dimer around, the quadratic can have numerical
    # problems.  But then we don't need to do quadratic anyway.  If [U] = Pt,
    # and some combination of Kd, Kf^2 and Dbar favored dimer, we would end up
    # with a more-than-physically-possible, conservation-of-mass-violating
    # concentration of dimer.  If this value (D_max) is still tiny, then there
    # isn't going to be any dimer around.  Get U, M, and return zeros for for
    # all other species.
    D_max = Pt*Kd*Kf*Kf*Dbar
    if _is_tiny(D_max,Pt):
        U = Pt/(1 + Kf)
        M = Kf*U
        return U, M, 0, 0, 0, 0, 0, 0

    # Quadratic coefficients in U
    a = 2*Kd*Kf*Kf*Dbar
    b = (1 + Kf)
    c = -Pt

    # Get positive root to find U
    U = (-b + np.sqrt(b**2 - 4*a*c))/(2*a)

    # If there is going to be U around, calculate M and D from U.
    if not _is_tiny(U,Pt):
        M = Kf*U
        D = Kd*(M**2)

    # If U is effectively zero, numerical problems likely.  Set U to 0 and do
    # quadratic in M.
    else:

        U = 0

        # Quadratic coefficients in M
        a = 2*Kd*Dbar
        b = 1
        c = -Pt

        # Get positive root to find M
        M = (-b + np.sqrt(b**2 - 4*a*c))/(2*a)

        # If there is going to be M around, calculate D from M.
        if not _is_tiny(M,Pt):
            D = Kd*(M**2)

        # If M is effectively zero, numerical problems likely.  Set M to zero
        # and solve for D directly
        else:
            M  = 0
            D = Pt/(2*Dbar)

    # Now get concentrations of dimer species given value for D we determined
    # above.
    Ds = Ks*D
    DsC1 = Ds*K1*C
    DsC2 = DsC1*K1*C
    DsC3 = DsC2*K2*C
    DsC4 = DsC3*K2*C

    # Return concentrations.  Note that the species DsC1 and DsC3 are possible
    # in two ways, so multiply by degeneracy of 2.
    return U, M, D, Ds, 2*DsC1, DsC2, 2*DsC3, DsC4


def _C_residual(params,Kf,Kd,Ks,K1,K2,Pt,Ct):
    """
    A residual function for finding the free value of calcium.  The
    value of free calcium is the fit parameter.  The difference between
    the calculated and total calcium and protein species' concentrations
    is the residual.
    """

    # Calculate the speccies concentration given our guess for [C]
    C = params[0]
    U, M, D, Ds, DsC1, DsC2, DsC3, DsC4 = _calc_prot_conc(Kf,Kd,Ks,K1,K2,Pt,C)

    # calculate the total concentrations in calcium and protein
    Ct_calc = C + DsC1 + 2*DsC2 + 3*DsC3 + 4*DsC4
    Pt_calc = U + M + 2*(D + Ds + DsC1 + DsC2 + DsC3 + DsC4)

    if Ct == 0:
        r_c = 0
    else:
        r_c = (Ct_calc - Ct)/Ct

    if Pt == 0:
        r_p = 0
    else:
        r_p = (Pt_calc - Pt)/Pt

    # Return difference between actual and total
    return np.array([r_c,r_p])

def _species_conc(Kf,Kd,Ks,K1,K2,
                 Pt,Ct,
                 convergence_cutoff=0.01,
                 guess_resolution=0.1,
                 verbose=False):
    """
    Get species of U, M, D, Ds, DsC1, DsC2, DsC3 and DsC4 given the equilibrium
    constants and the total protein and calcium concentrations. This function
    is private and assumes that Pt and Ct are both floats.  The public function
    does sanity checking and allows Pt and Ct to be arrays.

    Kf: folding equilibrium constant (Kf = M/U)
    Kd: dimerization equilibriucm constant (Kd=D/(M*M))
    Ks: apo to active dimer (Ks=Ds/D)
    K1: binding at stronger calcium binding site (K1=DsC1/(Ds*C)=DsC2/(DsC1*C))
    K2: binding at weaker calcium binding site (K2=DsC3/(DsC2*C)=DsC4/(DsC3*C))
    Pt: total protein monomer concentration
    Ct: total calcium concentration

    convergence_cutoff: percent difference between actual and calculated totals
    guess_resolution: if the first guess doesn't succeed, try another guess that
                      is guess-resolution away
    verbose: if True, record all all residuals and spit out if regression fails.

    returns: U, M, D, Ds, DsC1, DsC2, DsC3, DsC4
    """

    # If zero calcium total, we already know C
    if Ct == 0:
        C = 0.0
    else:

        best_residual = None
        best_C = None
        best_sum_residual = None

        guesses = np.arange(0,1+guess_resolution,guess_resolution)

        successful = False
        residuals = []
        for g in guesses:

            guess = Ct*g

            # Find value for C between 0 and Ct that finds species concentrations
            # that sum to Ct and Pt.
            fit = optimize.least_squares(_C_residual,
                                         [guess],
                                         bounds=(0,Ct),
                                         args=np.array((Kf,Kd,Ks,K1,K2,Pt,Ct)),
                                         tr_solver="exact",
                                         method="dogbox")

            # If we did not converge in the fit, try next guess
            if fit.success:
                C = fit.x[0]
            else:
                continue

            # get residuals
            residual = _C_residual([C],Kf,Kd,Ks,K1,K2,Pt,Ct)

            # Record residuals if this is verbose
            if verbose:
                residuals.append((g,C,residual))
                if best_sum_residual is None or np.sum(np.abs(residual)) < best_sum_residual:
                    best_residual = np.copy(residual)
                    best_sum_residual = np.sum(residual)
                    best_C = C

            # Have we converged?  If not, change the guess
            if np.abs(residual[0]) > convergence_cutoff or \
               np.abs(residual[1]) > convergence_cutoff:
                continue

            # If we get here, we've converged
            successful = True
            break

        if not successful:

            err = "\nCould not find solution within convergence cutoff.\n"
            if not verbose:
                err += "Run again, setting verbose = True to get information\n"
                err += "about the regression that is failing.\n"
            else:
                err += "Pt: {}\n".format(Pt)
                err += "Ct: {}\n".format(Pt)
                err += "Kf: {}\n".format(Kf)
                err += "Kd: {}\n".format(Kd)
                err += "Ks: {}\n".format(Ks)
                err += "K1: {}\n".format(K1)
                err += "K2: {}\n".format(K2)
                err += "\nAll residuals:\n\n"
                for r in residuals:
                    err += "{}\n".format(r)

                err += "\nBest residual:\n\n"
                err += "residual: {}\n".format(best_residual)
                err += "free [ca]: {}".format(best_C)

            raise RuntimeError(err)

    # Calculate concentration of each species given the equilibrium constants,
    # protein, and free calcium conc
    U, M, D, Ds, DsC1, DsC2, DsC3, DsC4 = _calc_prot_conc(Kf,Kd,Ks,K1,K2,Pt,C)

    # Make sure [U] is positive.
    if U < 0:
        err = "No positive value for [U] was found\n"
        raise RuntimeError(err)

    return U, M, D, Ds, DsC1, DsC2, DsC3, DsC4, C


def species_conc(Kf,Kd,Ks,K1,K2,
                 Pt,Ct,
                 convergence_cutoff=0.01,
                 guess_resolution=0.1,
                 verbose=False):
    """
    Get species of U, M, D, Ds, DsC1, DsC2, DsC3 and DsC4 given the equilibrium
    constants and the total protein and calcium concentrations.  This is the
    core, publically implementation of the model.  The equilibrium constants and
    concentrations may be floats or arrays, where all arrays have the same length.

    Kf: folding equilibrium constant (Kf = M/U)
    Kd: dimerization equilibriucm constant (Kd=D/(M*M))
    Ks: apo to active dimer (Ks=Ds/D)
    K1: binding at stronger calcium binding site (K1=DsC1/(Ds*C)=DsC2/(DsC1*C))
    K2: binding at weaker calcium binding site (K2=DsC3/(DsC2*C)=DsC4/(DsC3*C))
    Pt: total protein monomer concentration
    Ct: total calcium concentration

    convergence_cutoff: percent difference between actual and calculated totals
    guess_resolution: if the first guess doesn't succeed, try another guess that
                      is guess-resolution away
    verbose: if True, record all all residuals and spit out if regression fails.

    returns: U, M, D, Ds, DsC1, DsC2, DsC3, DsC4
    """

    mismatch_error = False
    values = [Kf,Kd,Ks,K1,K2,Pt,Ct]
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
        err = "Kf, Kd, Ks, K1, K2, Pt and Ct must either be float values or\n"
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
    Kf, Kd, Ks, K1, K2, Pt, Ct = values


    U, M, D, Ds, DsC1, DsC2, DsC3, DsC4, C = [], [], [], [], [], [], [], [], []
    for i in range(len(Pt)):

        concs = _species_conc(Kf[i],Kd[i],Ks[i],K1[i],K2[i],
                              Pt[i],Ct[i],
                              convergence_cutoff=convergence_cutoff,
                              guess_resolution=guess_resolution,
                              verbose=verbose)
        U.append(concs[0])
        M.append(concs[1])
        D.append(concs[2])
        Ds.append(concs[3])
        DsC1.append(concs[4])
        DsC2.append(concs[5])
        DsC3.append(concs[6])
        DsC4.append(concs[7])
        C.append(concs[8])

    U = np.array(U)
    M = np.array(M)
    D = np.array(D)
    Ds = np.array(Ds)
    DsC1 = np.array(DsC1)
    DsC2 = np.array(DsC2)
    DsC3 = np.array(DsC3)
    DsC4 = np.array(DsC4)
    C = np.array(C)

    return U, M, D, Ds, DsC1, DsC2, DsC3, DsC4, C
