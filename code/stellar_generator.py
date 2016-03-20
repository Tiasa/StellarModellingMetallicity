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
readable_strings = ("Radius", "Density", "Temperature", "Luminosity", "Mass", "Opt Depth")

class Star():

    def __init__(self, temp_c, density_c, composition):
        self.composition = composition
        self.temp_c = temp_c
        self.density_c = density_c
        self.__solved = False

    def diP_diT(self, ss, r):
        """
        The partial pressure gradient with respect to temperature
        """
        ideal_gas = ss[density] * k / (self.composition.mu * m_p)
        photon_gas = 4/3 * a * ss[temp]**3
        return ideal_gas + photon_gas

    def diP_dirho(self, ss, r):
        """
        The partial pressure gradient with respect to pressure
        """
        ideal_gas = k * ss[temp] / (self.composition.mu * m_p)
        nonrel_degenerate = nonrelgenpress * ss[density]**(2/3)
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
        kappa = 1/((1/kH)+(1/max(kff,kes)))
        return kappa

    def pressure(self, ss, r):
        """
        Pressure at a given state in the star due to three competing effects
        """
        nonrel_degenerate = nonrelgenpress * ss[density]**(5/3)
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
        radiative = 3 * self.opacity(ss, r) * ss[density] * ss[lumin] / ( 16 * pi * a * c * ss[temp]**3 * r**2)
        convective = ( 1 - 1 / gamma) * ( ss[temp] / self.pressure(ss, r) ) * ( G * ss[mass] * ss[density] ) / ( r**2 )
        return -min(radiative, convective)

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

    def solve(self):
        if self.__solved:
            raise Exception("Star already solved.")

        r_0 = 1

        ic = np.empty(ss_size)

        ic[density] = self.density_c
        ic[temp] = self.temp_c
        ic[mass] = 4 * pi / 3 * r_0**3 * self.density_c
        ic[lumin] = ic[mass] * self.density_c * self.energy_prod(ic, r_0)
        ic[opt_depth] = 0 # ???

        system_DE = [self.drho_dr, self.dT_dr, self.dL_dr, self.dM_dr, self.dtau_dr] # Needs to match stellar state enum order above
        R = 7e9
        tolerance = 100
        max_step = 1e10
        min_step = 1e-10

        r, ss = rkf(system_DE, r_0, R, ic, tolerance, max_step, min_step)

        self.r_profile = r
        self.ss_profile = ss
        self.__solved = True

    def plot(self):
        if self.__solved is False:
            self.solve()

        for i in range(self.ss_profile.shape[0]):
            plt.figure()
            plt.title("{0} vs. Radius".format(readable_strings[i+1]))
            plt.xlabel("Radius")
            plt.ylabel(readable_strings[i+1])
            plt.plot(self.r_profile, self.ss_profile[i, :])
            # plt.plot(self.r_profile[:250], self.ss_profile[i, :][:250])
            plt.show()

    def plot_step_sizes(self):
        if self.__solved is False:
            self.solve()

        steps = self.r_profile[1:] - self.r_profile[:-1]
        plt.figure()
        plt.plot(range(len(steps)), steps)
        plt.show()

    def log(self, a=0, b=0):
        if self.__solved is False:
            self.solve()
        size = self.ss_profile.shape[1]
        print("{0:20} | {1:20} | {2:20} | {3:20} | {4:20} | {5:20}".format(*readable_strings))
        log_format = "{0:20.10E} | {1:20.10E} | {2:20.10E} | {3:20.10E} | {4:20.10E} | {5:20.10E}"
        if a>0:
            print("Printing first {0} data entires.".format(a))
            for i in np.arange(0, min(size, a), 1):
                print(log_format.format(self.r_profile[i], *self.ss_profile[:, i]))
        if b>0:
            print("Printing last {0} data entires.".format(b))
            for i in np.arange(min(size, size-b), size, 1):
                print(log_format.format(self.r_profile[i], *self.ss_profile[:, i]))


# test_star = Star(temp_c = 1.5e7, density_c=1.6e5, composition=Composition.fromXY(0.69, 0.29))
test_star = Star(temp_c = 1.5e7, density_c=1.6e5, composition=Composition.fromXY(0.69, 0.29))

test_star.solve()
test_star.plot()
test_star.plot_step_sizes()
test_star.log(b=20)
