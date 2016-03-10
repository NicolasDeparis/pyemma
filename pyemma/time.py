# coding: utf-8

import numpy as np
from scipy import integrate

def a2t_quad(az, o_m, H0):

    o_v=1.-o_m

    H0 = H0*1e3/1e6/3.08567758e16*(365*24*3600) # Hubble constant in SI unit

    def f(a):
        return a/ np.sqrt(o_m*a + o_v*a**4 )

    return 1./H0 * integrate.quad(f,0,az)[0]
