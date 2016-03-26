from __future__ import division, print_function
import numpy as np

## All the contants
G = 6.674e-11 # Gravitational Constant
c = 2.998e8   # Speed of light
h = 6.626e-34 # Planck's constant
hbar = h / (2 * np.pi)
k = 1.381e-23 # Boltzmann Constant
sigma = 5.670e-8 # Stefan-Boltzmann Constant
m_e = 9.109e-31 # Electron mass
M_s = 1.989e30 # Mass of Sun
R_s = 6.963e8 # Radius of Sun
L_s = 3.839e26 # Luminosity of Sun
m_p = 1.6726219e-27 # Mass of proton
a = 4 * sigma / c
pi = np.pi
ep_pp_coeff = 1.07e-7 * 1e-5 * (1e-6)**4
ep_cno_coeff = 8.24e-26 * 1e-5 * (1e-6)**19.9
nonrelgenpress = (3 * pi)**(2/3)/5 * hbar**2/m_e * m_p**(-5/3)
mach_ep = np.finfo(np.float64).eps
tiny_float = 1e-50
gamma = 5/3 # ideal gas constant
