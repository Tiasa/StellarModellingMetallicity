from __future__ import division, print_function

from constants import *
from energy_production import *
from stellar_structure import *
from composition import Composition
from rkf import rkf
import matplotlib.pyplot as plt

# --- Stellar State Enum ---
ss_size = 5
# --------------------------
density = 0
temp = 1
lumin = 2
mass = 3
opt_depth = 4

class Star():

    def __init__(self, temp_c, density_c, composition):
        self.composition = composition
        self.temp_c = temp_c
        self.density_c = density_c
        self.__solved = False

    def solve(self):
        if self.__solved:
            raise Exception("Star already solved.")

        r_0 = mach_ep

        ic = np.empty(ss_size)

        ic[density] = self.density_c
        ic[temperature] = self.temp_c
        ic[mass] = 4 * pi / 3 * r_0**3 * self.density_c
        ic[lumin] = ic[mass] * self.density_c * self.energy_prod(ic, r_0)
        ic[opt_depth] = 0 # ???

        system_DE = [drho_dr, dT_dr, dL_dr, dM_dr, dtau_dr] # Needs to match stellar state enum order above
        R = 10000
        tolerance = 1e-6
        max_step = 100
        min_step = 1

        r, ss = rkf(system_DE, r_0, R, ic, tolerance, max_step, min_step)

        self.r_profile = r
        self.ss_profile = ss
        self.__solved = True

    def diP_diT(self, ss, r):
        """
        The partial pressure gradient with respect to temperature
        """
        ideal_gas = ss[density] * k / (self.composition.mu * m_p)
        photon_gas = 4/3 * a * ss[temp]**3
        return ideal_gas + radiative

    def diP_dirho(self, ss, r):
        """
        The partial pressure gradient with respect to pressure
        """
        ideal_gas = k * ss[temp] / (self.composition.mu * m_p)
        nonrel_degenerate = nonrelgenpress * (ss[density]/m_p)**(2/3)
        return ideal_gas + nonrel_degenerate

    def opacity(self, ss, r):
        """
        Opacity of the star as it depends on composition, density and temperature
        """
        ## Electron Scattering Opacity
        kes = 0.2 * (1+self.composition.X)
        ## Free-free Opacity
        kff = 1e24 * (self.composition.Z+0.0001) * (ss[density]**0.7) * (ss[temp]**(-3.5))
        ## Hydrogen opacity
        kH = 2.5e-32 * (self.composition.Z/0.02) * (ss[density]**0.5) * (ss[temp]**9)
        ## Kappa
        kappa = (1/kH)+(1/np.max(kff,kes))
        return kappa

    def pressure(self, ss, r):
        """
        Pressure at a given state in the star due to three competing effects
        """
        nonrel_degenerate = nonrelgenpress * (ss[density]/m_p)**(5/3)
        ideal_gas = k * ss[density] * ss[temp] / (self.composition.mu * m_p)
        photon_gas = 4/3 * a * ss[temp]**4
        return nonrel_degenerate + ideal_gas + photon_gas

    def energy_prod(self, ss, r):
        """Total energy production for a given X, density and temperature"""
        X = self.composition.X
        ep_pp = ep_pp_coeff * ss[density] * X**2 * ss[temp]**4
        ep_cno = ep_cno_coeff * ss[density] * X * (X * 0.03) * ss[temp]**19.9 # X_CNO = X * 0.03?
        return ep_pp + ep_cno

    # Stellar Structure
    def dT_dr(self, ss, r):
        """
        Temperature gradient with respect to radius
        """
        dT_dr_1 = 3 * self.opacity(ss, r) * ss[density] * ss[lumin] / ( 16 * pi * a * c * ss[temp]**3 * r**2)
        dT_dr_2 = ( 1 - 1 / gamma) * ( ss[temp] / self.pressure(ss, r) ) * ( G * ss[mass] * ss[density] ) / ( r**2 )
        return -min(dT_dr_1, dT_dr_2)

    # Stellar Structure
    def drho_dr(self, ss, r):
        """
        Density gradient with respect to radius
        """
        return - ((G * ss[mass] * ss[density])/(r**2) + self.diP_diT(ss,r) * self.dT_dr(ss, r))/(self.diP_dirho(ss,r))

    # Stellar Structure
    def dM_dr(self, ss, r):
        """
        Mass gradient with respect to radius
        """
        return 4 * pi * r**2 * ss[density]

    # Stellar Structure
    def dL_dr(self, ss, r):
        """
        Luminosity gradient with respect to radius
        """
        return 4 * pi * r**2 * ss[density] * self.energy_prod(ss, r)

    # Stellar Structure
    def dtau_dr(self, ss, r):
        """
        Optical depth gradient
        """
        return self.opacity(ss, r) * ss[density]


test_star = Star(temp_c = 1.5e7, density_c=1.6e5, composition=Composition.fromXY(0.69, 0.29))

test_star.solve()

# def stellar_solver(T_c, rho_c, X, Y):

#     buffer = 1e5



#     r          = np.empty(buffer)
#     r[[0,1]]   = r_i
#     rho        = np.empty(buffer)
#     rho[[0,1]] = rho_i
#     T          = np.empty(buffer)
#     T[[0,1]]   = T_i
#     M          = np.empty(buffer)
#     M[[0,1]]   = M_i
#     L          = np.empty(buffer)
#     L[[0,1]]   = L_i

#     i = 0
#     initial_step = 1e5
#     step = initial_step
#     while i < 100:
#         print(rkf(func, 3, i, 1))
#         i += 1
#     pass


# def fa(x,t):
#     return x[1]

# def fb(x,t):
#     return -x[0]

# def fc(x,t):
#     return x[3] + x[0]

# def fd(x,t):
#     return -x[2] - x[1] * x[0]

# T, X = rkf(np.array([fa,fb, fc, fd]), 0, 2 * pi, np.array([1, 0, -1, 0.3]), )
# # print(T.shape, X.shape, X[1, :].shape)

# plt.figure()
# for i in range(X.shape[0]):
#     plt.plot(T, X[i, :])
# plt.show()


# # Tests
# # stellar_solver(
# # )