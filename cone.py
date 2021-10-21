import numpy as np


def compute_reflected_light(Iin, a_cone, p_conc, a_scatter=1):
    """
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
    """
    absorption_fraction = a_cone * p_conc
    attenuation = 1 - absorption_fraction
    return a_scatter * Iin * attenuation


def dpdt(Iin, Qe, p_conc, t0=120):
    """
    Compute the derivative of the pigment concentration at a given state.

    Parameters
    ----------
    Iin : scalar (1,)
    Input image intensity at a given time
    Qe : scalar (1,)
    Time constant in units of td-sec
    p_conc : scalar (1,)
    Concentration of photopigment at a given time
    """
    return (-Iin * p_conc / Qe) + ((1 - p_conc) / t0)


def cone_L(lmbda):
    """
    Interpolate the response of the L cone for a given wavelength lmbda.
    """
    xp = [400, 475, 500, 525, 550, 575, 600, 625, 650, 675, 700]
    fp = [0, 0.1, 0.3, 0.7, 0.95, 0.999, 0.9, 0.553, 0.18, 0.075, 0]
    alpha_cone = np.interp(lmbda, xp, fp)
    return alpha_cone


def cone_M(lmbda):
    """
    Interpolate the response of the M cone for a given wavelength lmbda.
    """
    xp = [400, 450, 500, 525, 550, 575, 600, 625, 650, 675, 700]
    fp = [0, 0.1, 0.4, 0.9, 1, 0.8, 0.35, 0.1, 0.05, 0.001, 0]
    alpha_cone = np.interp(lmbda, xp, fp)
    return alpha_cone


def cone_S(lmbda):
    """
    Interpolate the response of the S cone for a given wavelength lmbda.
    """
    xp = [400, 425, 447, 475, 500, 512.5, 550]
    fp = [0.05, 0.6, 0.998, 0.45, 0.1, 0.05, 0]
    alpha_cone = np.interp(lmbda, xp, fp)
    return alpha_cone
