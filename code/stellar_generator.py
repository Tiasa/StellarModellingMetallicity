from __future__ import division, print_function

from constants import *
from energy_production import *

def stellar_solver(T_c, rho_c):
    pass

def diP_diT(rho, mu, T):
    """
    The partial pressure gradient with respect to temperature
    """
    ideal_gas = rho * k / (mu * m_p)
    photon_gas = 4/3 * a * T**3
    return ideal_gas + radiative

def diP_dirho(rho, mu, T):
    """
    The partial pressure gradient with respect to pressure
    """
    ideal_gas = k * T / (mu * m_p)
    nonrel_degeneate = nonrelgenpress * rho**(2/3)
    return ideal_gas + nonrel_degeneate

def drho_dr(M, rho, r, diP_diT, dT_dr, diP_dirho):
    """
    Density gradient
    """
    return - ((G * M * rho)/(r**2) + diP_diT * dT_dr)/(diP_dirho)

def dM_dr(r,rho):
    """
    Mass Gradient
    """
    return 4 * pi * r**2 * rho

