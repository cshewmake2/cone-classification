import numpy as np
import cone


def estimate_a_scatter(duration, a_scatter, a_cone, Iin, mu_m, sigma_m):
    '''
    Estimate a_intrins, intrinsic refectance assume imaging with a wavelength 
    outside of visual spectrum (no bleaching)

    Parameters
    ----------
    duration : scalar (1,), units = sec
        Time to infer intrinsic reflectance over
    a_scatter : scalar (1,), default = 1
        Fraction of light remaining after scattering, unique to each cone
    a_cone : scalar (1,) in range [0,1]
        Cone-specific absorption fraction at the stimulus frequency
    Iin : scalar (1,), units = td-sec
        Applied retinal illuminance by AOSLO
    mu_m : scalar (1,)
        Mean value of measurement noise
    sigma_m : scalar (1,)
        Std of measurement noise

    Return(s)
    ---------
    a_scatter_est : scalar (1,)
    Estimated a_scatter
    '''
    update_rate = 1/30    # [s] update rate of Iin(t) (AOSLO fps)
    
    # bleach each cones
    pt = 0 

    t = np.arange(0, duration, update_rate)
    Iout = []
    

    for i, ti in enumerate(t):
        
        # noise/scatter model
        nfl_noise = np.random.normal(loc=mu_m,scale=sigma_m)
        Iout_ti = cone.compute_reflected_light(Iin,            \
                                          a_cone,              \
                                          pt,                  \
                                          a_scatter=a_scatter, \
                                         ) + nfl_noise
        Iout.append(Iout_ti)
    return np.mean((np.asarray(Iout)-mu_m)/Iin)


def simulate(**kwargs):
    '''
    Simulate optical feedback method

    Parameters
    ---------- 
    account_for_a_scatter : boolean, default=False
        Use method to estimate scatter and correct in Iin(t+1) calculation
    update_rate : scalar (1,) default = 0.1 (10Hz)
        Update rate of Iin(t) = Iout(t-1)
    dt : scalar (1,) default = 0.001 (1ms)
        Step size of cone bleaching dynamics model
    n_seconds : scalar (1,) default = 6 (sec)
        Duration to stimulate experiment
    scatter_est_duration : scalar (1,) default = 5 (sec)
        Time to estimate a_scatter parameter
    p0 : scalar (1,) default = 1. (dark adapted)
        Initial value of cone bleach
    Qe : scalar (1,) default = 3e6 (td-sec)
        Time constant of bleaching dynamics in dark-bright conditions
    mu_m : scalar (1,) default = 250
        Mean value of measurement noise
    sigma_m : scalar (1,) default = 500
        Std of measurement noise
    Iin_0 : scalar (1,) default = 4e4 (td-sec)
        Initial input light intensity
    a_scatter_bounds : scalar tuple (2,) default = (1,1)
        defines low & high values of rand. uniform dist. a_scatter is drawn from
        
    Return(s)
    ---------
    Iin : scalar (n_seconds/dt+1,)
        Input light over time
    Iout : scalar (n_seconds/dt,)
        Output light over time
    pt : scalar (n_seconds/dt+1,)
        Bleaching dynamics throughout experiment
    t : scalar (n_seconds/dt,)
        time
    a_scatter_residual : scalar (1,)
        residiual between a_scatter and estimate
    '''
    # set defaults 
    a_cone=kwargs.pop('a_cone',0.1)
    account_for_a_scatter=kwargs.pop('account_for_a_scatter',False)
    update_rate=kwargs.pop('update_rate',0.1)
    dt=kwargs.pop('dt',0.001)
    n_seconds=kwargs.pop('n_seconds',6)
    scatter_est_duration=kwargs.pop('scatter_est_duration',5)
    p0=kwargs.pop('p0',1.)
    Qe=kwargs.pop('Qe',3e6)
    mu_m=kwargs.pop('mu_m',250)
    sigma_m=kwargs.pop('sigma_m',500)
    Iin_0=kwargs.pop('Iin_0',4e4)
    a_scatter_bounds=kwargs.pop('a_scatter_bounds',(1,1))
    
    
    a_scatter = np.random.uniform(low=a_scatter_bounds[0], high=a_scatter_bounds[1])
  
    t = np.arange(0, n_seconds, dt)
    update_t = np.arange(update_rate, n_seconds, update_rate)
    
    Iin = [Iin_0]
    Iout = []

    pt = [p0]
    p_ti = p0

    # estimate a_scatter
    a_scatter_residual = 9e9
    if account_for_a_scatter:
        a_scatter_est = estimate_a_scatter(duration=scatter_est_duration, \
                                             a_scatter=a_scatter,         \
                                             a_cone=a_cone,               \
                                             Iin=Iin_0,                   \
                                             mu_m=mu_m,                   \
                                             sigma_m=sigma_m              \
                                          )
    
        a_scatter_residual = np.abs(a_scatter_est-a_scatter)
        # print('abs difference between known and estimated intrinsic reflectance: %.4e' % (a_intrins_residual))

    updatei = 0
    for i, ti in enumerate(t):
        # update photopigment concentration
        p_grad = cone.dpdt(Iin[i], Qe, pt[i])
        p_ti += p_grad * dt
        pt.append(p_ti)

        # noise/scatter model
        nfl_noise = np.random.normal(loc=mu_m,scale=sigma_m)

        Iout_ti = cone.compute_reflected_light(Iin[i],         \
                                          a_cone,              \
                                          p_ti,                \
                                          a_scatter=a_scatter, \
                                         ) + nfl_noise
        Iout.append(Iout_ti)

        # check if its time to update Iin
        if updatei < update_t.size and 0 <= (ti-update_t[updatei]) < dt:
            # update input image from output image
            if account_for_a_scatter:
                Iout_corrected = (Iout[i]-mu_m)/a_scatter_est
                Iin.append(np.clip(Iout_corrected, 0, 40000))
            else:
                Iin.append(np.clip(Iout[i], 0, 40000))
            updatei = updatei+1
        else:
            # duplicate the previous Iin
            Iin.append(Iin[i])

    return np.asarray(Iin), np.asarray(Iout), np.asarray(pt), t, a_scatter_residual
