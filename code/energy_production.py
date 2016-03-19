from __future__ import division, print_function
from constants import *

"""
Energy production calculations
"""

def ep_pp(rho, X, T):
    """Proton-proton chain energy production"""
    return ep_pp_coeff * rho * X**2 * T**4

def ep_cno(rho, X, T):
    """CNO cycle energy production"""
    return ep_cno_coeff * rho * X * X_CNO(X) * T**19.9

def ep(rho, X, T):
    """Total energy production for a given X, density and temperature"""
    return ep_pp(rho, X, T) + ep_cno(rho, X, T)

def X_cno(X):
    return X * 0.03