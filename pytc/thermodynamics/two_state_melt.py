

def species_conc(Ku,Pt):
    """
    Two-state unfolding.

    Ku: equilibrium constant for unfolding
    Pt: total protein concentration

    Returns:
    F: concentration of folded
    U: concentration of unfolded
    """

    U = 1/(1+Ku)*Pt
    F = Pt - U

    return F, U
