from __future__ import division
from __future__ import print_function

_ep_pp_coeff = 1.07e-7 * 1e-5 * (1e-6)**4
_ep_cno_coeff = 8.24e-26 * 1e-5 * (1e-6)**19.9

def ep_pp(rho, X, T):
    """Proton-proton chain energy production"""
    return _ep_pp_coeff * rho * X**2 * T**4

def ep_cno(rho, X, T):
    """CNO cycle energy production"""
    return _ep_cno_coeff * rho * X * X_CNO(X) * T**19.9

def ep(rho, X, T):
    """Total energy production for a given X, density and temperature"""
    return ep_pp(rho, X, T) + ep_cno(rho, X, T)

def X_cno(X):
    return X * 0.03