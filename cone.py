import numpy as np
from scipy import interpolate


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


class cone():
    def __init__(self,fname):
        '''
        Class for cone types
        Parameters
        ----------
        fname : string 
            path to csv file that contains cone activation as function of lmbda
        '''
        data = np.genfromtxt(fname, delimiter=',')
        data[np.isnan(data)] = 0.
        self.lmbda = data[:,0]
        self.L = data[:,1]
        self.M = data[:,2]
        self.S = data[:,3]


    def cone_L(self,querylmbda):
        """
        Interpolate the response of the L cone for a given wavelength lmbda.
        """
        alpha_cone = np.interp(querylmbda, self.lmbda, self.L)
        alpha_cone = np.clip(alpha_cone, a_min=0.0, a_max=1)
        return alpha_cone


    def cone_M(self,querylmbda):
        """
        Interpolate the response of the M cone for a given wavelength lmbda.
        """
        alpha_cone = np.interp(querylmbda, self.lmbda, self.M)
        alpha_cone = np.clip(alpha_cone, a_min=0.0, a_max=1)
        return alpha_cone


    def cone_S(self,querylmbda):
        """
        Interpolate the response of the S cone for a given wavelength lmbda.
        """
        alpha_cone = np.interp(querylmbda, self.lmbda, self.S)
        alpha_cone = np.clip(alpha_cone, a_min=0.0, a_max=1)
        return alpha_cone
