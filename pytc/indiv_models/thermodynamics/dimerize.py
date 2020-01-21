
def species_conc(Kdimer,Pt):
    """
    """

    a = 2*Kdimer
    b = 1
    c = -Pt

    M = (-b + np.sqrt(1 + 8*Kdimer*Pt))/(4*Kdimer)
    D = (Pt - M)/2

    return M, D

def signal(param,df,prot_conc_key="Pt"):

    Kdimer = param[0]
    signal_M = param[1]
    signal_D = param[2]

    M, D = species_conc(Kdimer,df[prot_conc_key])

    return M*signal_M + D*signal_D
