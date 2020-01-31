
def species_conc(Kdimer,Pt):
    """
    """

    a = 2*Kdimer
    b = 1
    c = -Pt

    M = (-b + np.sqrt(1 + 8*Kdimer*Pt))/(4*Kdimer)
    D = (Pt - M)/2

    return M, D
