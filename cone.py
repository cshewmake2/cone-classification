import numpy as np

def compute_reflected_light(Iin, a_cone, p_conc, a_scatter = 1):
    '''
    Compute the light reflected for a given cone type & pigment concentration.

    Parameters
    ----------
    Iin : scalar (1,) 
    Input image intensity at a given time
    a_cone : scalar (1,) in range [0,1]
    Cone-specific absorption fraction at the stimulus frequency
    p_conc : scalar (1,)
    Concentration of photopigment at a given time
    a_scatter : scalar (1,), default = 1
    Fraction of light remaining after scattering
    a_intrins : scalar (1,), default = 1
    Intrinsic reflectance of cone
    '''
    absorption_fraction = a_cone * p_conc
    attenuation = 1 - absorption_fraction
    return a_scatter * Iin * attenuation

def dpdt(Iin, Qe, p_conc, t0=120): 
    '''
    Compute the derivative of the pigment concentration at a given state.

    Parameters
    ----------
    Iin : scalar (1,) 
    Input image intensity at a given time
    Qe : scalar (1,)
    Time constant in units of td-sec
    p_conc : scalar (1,)
    Concentration of photopigment at a given time
    '''
    return (-Iin * p_conc / Qe) + ((1 - p_conc)/t0)
