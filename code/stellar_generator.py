from __future__ import division, print_function

from constants import *
from energy_production import *

def stellar_solver(T_c, rho_c, X, Y):

    r_0 = mach_ep

    r_i = [0, r_0]
    rho_i = [rho_c, rho_c]
    T_i = [T_c, T_c]
    M_i = [0, 4 * pi / 3 * r_0**3 * rho_c]
    L_i = [0, M_i[1] * rho_c * ep(rho_c, X, T_c) ]

    r = np.array(r_i)
    M = np.array(M_i)
    rho = np.array(rho_i)
    T = np.array(T_i)
    M = np.array(M_i)


    while i < 100:

        i++
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

def dT_dr(kappa, rho, L, T, P, M, r):
    """
    Temperature gradient
    """
    dT1 = abs( 3 * kappa * rho * L / ( 16 * pi * a * c * T**3 * r**2) )
    dT2 = abs( ( 1 - 1 / gamma) * ( T / P ) * ( G * M * rho ) / ( r**2 ) )
    return min(dT1, dT2)

def dL_dr(r, rho, ep):
    """
    Luminosity gradient
    """
    return 4 * pi * r**2 * rho * ep

def dtau_dr(kappa, rho):
    """
    Optical depth gradient
    """
    return kappa * rho
