

def simulate_expt(Kd,Pt,Xt):

    A, B, AB = single_site_binding.species_conc(Kd,Pt,Xt)

    return AB/At


def yo(param,df,observe_P=True):
    """
    """

    Kd = param[0]

    P, X, PX = single_site_binding.species_conc(Kd,df.Pt,df.Xt)

    if observe_P:
        return PX/df.Pt
    else:
        return PX/df.Xt

def yo2(param,df,titrant_key="osm",prot_conc_key=None):

    dG0 = param[0]
    m = param[1]

    Ku = dG0 + m*df[titrant_key]
    if prot_conc_key is None:
        prot_conc = df[prot_conc_key]
    else:
        prot_conc = np.ones(len(Ku))


    F, U = two_state_melt.species_conc(Ku,prot_conc)

    return F/(F + U)


def yo3(param,df):

    Kdimer = param[0]

    
