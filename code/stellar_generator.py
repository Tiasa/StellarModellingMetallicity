from __future__ import division, print_function

from constants import *
from composition import Composition
from rkf import rkf
from adaptive_bisection import adaptive_bisection

# --- Stellar State Enum ---
ss_size = 4
# --------------------------
density = 0
temp = 1
lumin = 2
mass = 3
# opt_depth = 4
readable_strings = ("Radius", "Density", "Temperature", "Luminosity", "Mass") #, "Opt Depth")

class Star():

    def __init__(self, temp_c, composition):
        self.composition = composition
        self.temp_c = temp_c
        self.delta_tau_thres = 1e-20
        self.stellar_structure_eqns = [self.drho_dr, self.dT_dr, self.dL_dr, self.dM_dr] # Needs to match stellar state enum order above
        self.stellar_structure_size = ss_size
        self.is_solved = False

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
        nonrel_degenerate = (5/3) * nonrelgenpress * ss[density]**(2/3)
        return ideal_gas + nonrel_degenerate

    def dP_dr(self, ss, r):
        return - (G * ss[mass] * ss[density])/(r**2)

    def partial_opacity(self, ss, r):
        # Electron Scattering Opacity
        kes = 0.02 * (1+self.composition.X)
        # Free-free Opacity
        kff = 1e24 * (self.composition.Z+0.0001) * ((ss[density] * 1e-3)**0.7) * (ss[temp]**(-3.5))
        # Hydrogen opacity
        kH = 2.5e-32 * ((self.composition.Z+tiny_float)/0.02) * ((ss[density] * 1e-3)**0.5) * (ss[temp]**9)
        # Avoid division by zero (kes, kff cant be zero) added tiny_float

        return (kes, kff, kH)

    def opacity(self, ss, r):
        """
        Opacity of the star as it depends on composition, density and temperature
        """
        (kes, kff, kH) = self.partial_opacity(ss,r)

        kappa = 1/((1/kH)+(1/max(kff,kes)))
        return kappa

    def partial_pressure(self, ss, r):
        """
        Each of the three pressure sources at different points in th star
        """
        nonrel_degenerate = nonrelgenpress * ss[density]**(5/3)
        ideal_gas = k * ss[density] * ss[temp] / (self.composition.mu * m_p)
        photon_gas = 1/3 * a * ss[temp]**4
        return (nonrel_degenerate, ideal_gas, photon_gas)

    def pressure(self, ss, r):
        """
        Pressure at a given state in the star due to three competing effects
        """
        return sum(self.partial_pressure(ss,r))

    def partial_energy_prod(self, ss, r):
        X = self.composition.X
        ep_pp = ep_pp_coeff * ss[density] * X**2 * ss[temp]**4
        ep_cno = ep_cno_coeff * ss[density] * X * (X * 0.03) * ss[temp]**19.9 # X_CNO = X * 0.03?
        return (ep_pp, ep_cno)

    def energy_prod(self, ss, r):
        """
        Total energy production for a given X, density and temperature
        """
        return sum(self.partial_energy_prod(ss,r))

    def delta_tau(self, ss, r):
        """
        Delta tau stopping condition used to determine when r is well outside the star
        """
        return (self.opacity(ss, r) * ss[density]**2) / abs(self.drho_dr(ss, r))

    def partial_dT_dr(self, ss, r):
        """
        The two methods of energy transport
        """
        radiative = (3 * self.opacity(ss, r) * ss[density] * ss[lumin]) / ( 16 * pi * a * c * ss[temp]**3 * r**2)
        convective = ( 1 - 1 / gamma) * ( ss[temp] / self.pressure(ss, r) ) * ( G * ss[mass] * ss[density] ) / ( r**2 )

        return (radiative, convective)

    # Stellar Structure
    def dT_dr(self, ss, r):
        """
        Temperature gradient with respect to radius
        """
        return -min(*self.partial_dT_dr(ss, r))

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

    def partial_dL_dr(self, ss, r):
        """
        Luminosity gradient with respect to radius
        """
        partial_energy_production = self.partial_energy_prod(ss,r)
        mass_grad = self.dM_dr(ss,r)
        dL_dr_pp = mass_grad * partial_energy_production[0]
        dL_dr_cno = mass_grad * partial_energy_production[1]
        return (dL_dr_pp, dL_dr_cno)

    # Stellar Structure
    def dL_dr(self, ss, r):
    	"""
    	Luminosity gradient with respect to radius
    	"""
    	return sum(self.partial_dL_dr(ss, r))

    # Stellar Structure
    def dtau_dr(self, ss, r):
        """
        Optical depth gradient
        """
        return - self.opacity(ss, r) * ss[density]

    def solve(self):
        if self.is_solved:
            return

        # (i_surf, ss, r, delta_tau_surf) = self.solve_density_c(1.5e3)
        # raise Exception("fds")
        # density_c = adaptive_bisection(self.solve_density_c_error, 1e-1, 1e14, 0.01)
        # density_c = adaptive_bisection(self.solve_density_c_error, 1e4, 1e9, 0.01)
        density_c, tol = adaptive_bisection(self.solve_density_c_error, 1, 1e8)
        # density_c = adaptive_bisection(self.solve_density_c_error, 0.03, 500, 1)

        # print("---- Solving Star With Correct Central Density ---")
        (i_surf, ss, r, delta_tau_surf) = self.solve_density_c(density_c, tol)
        # print("--------------------- Solved ---------------------")

        self.i_surf = i_surf
        self.ss_profile = ss
        self.r_profile = r
        self.ss_surf = ss[:,i_surf]
        self.r_surf = r[i_surf]
        self.density_c = density_c
        self.delta_tau_surf = delta_tau_surf
        self.lumin_surf_bb, self.lumin_surf_rkf = self.relative_surface_lumin(i_surf, ss, r)
        self.temp_surf = (self.lumin_surf_rkf / (4 * pi * sigma * self.r_surf**2))**(1/4)
        self.lumin_surf = self.ss_surf[lumin]
        self.mass_surf = self.ss_surf[mass]
        self.density_surf = self.ss_surf[density]
        self.data_size = len(r)

        self.is_solved = True

    def relative_surface_lumin(self, i, ss, r):
        lumin_surf_bb = 4 * pi * sigma * r[i]**2 * ss[temp, i]**4
        lumin_surf_rkf = ss[lumin, i]

        return lumin_surf_bb, lumin_surf_rkf

    def solve_density_c_error(self, density_c, tol):
        lumin_surf_bb, lumin_surf_rkf = self.relative_surface_lumin(*self.solve_density_c(density_c, tol)[0:3])

        error = (lumin_surf_rkf - lumin_surf_bb)/np.sqrt(lumin_surf_rkf * lumin_surf_bb)
        return error

    def system_DE(self, ss, r):
        if np.min(ss) < 0:
            raise ValueError("Stellar state values are negative.")
        # elif np.isinf(self.opacity(ss, r)) or np.isnan(self.opacity(ss, r)):
        #     raise ValueError("Opacity too big.")
        else:
            return np.array([f(ss,r) for f in self.stellar_structure_eqns])

    def solve_density_c(self, density_c, tol):
        if self.is_solved:
            return

        temp_c = self.temp_c

        r_0 = 1 # Doesn't matter really as long as this is less than scale height

        ic = np.empty(self.stellar_structure_size)

        ic[density] = density_c
        ic[temp] = temp_c
        ic[mass] = 4 * pi / 3 * r_0**3 * density_c
        ic[lumin] = ic[mass] * self.energy_prod(ic, r_0)

        # print(self.delta_tau_thres)
        tau_surf = 2/3

        r, ss = rkf(self.system_DE, r_0, ic, tol, self.stop_condition)

        surfaced = False

        for i in reversed(xrange(1, r.shape[0]-1)):
            i_delta_tau_inner = self.delta_tau(ss[:, i-1], r[i-1])
            i_delta_tau_outer = self.delta_tau(ss[:, i], r[i])
            if (i_delta_tau_inner > tau_surf):
                if abs(i_delta_tau_inner - tau_surf) > abs(i_delta_tau_outer - tau_surf):
                    i_surf = i # Take outer
                else:
                    i_surf = i-1 # Take inner

                i_surf = i_surf
                ss_surf = ss[:, i_surf]
                r_surf = r[i_surf]
                delta_tau_surf = self.delta_tau(ss[:, i_surf], r[i_surf])
                surfaced = True
                break

        assert surfaced, "Surface wasn't reached"

        return (i_surf, ss, r, delta_tau_surf)

    def stop_condition(self, i, ss, r):
        delta_tau = self.delta_tau(ss, r)
        mass_limit = 1e3 * M_s
        # if ss[mass] > mass_limit:
            # print(delta_tau)
            # print("Terminating star due to mass limit.")
            # return True
        if delta_tau < self.delta_tau_thres:
            # print(i,ss,r)
            # print(delta_tau)
            # print("Terminating star due to optical depth.")
            return True
        return False

    def log_solved_properties(self):
        print("------- Solved Variables -------")
        log_format = "{0:40} {1:10.10E}"
        print(log_format.format("Delta tau surface (2/3 = 6.7E-1)", self.delta_tau_surf))
        print(log_format.format("Central density", self.density_c))
        print(log_format.format("Luminosity Blackbody", self.lumin_surf_bb))
        print(log_format.format("Luminosity RKF", self.lumin_surf_rkf))
        print(log_format.format("Inferred Surface Temperature", self.temp_surf))
        print(log_format.format("Data Size", self.data_size))

        self.log_ss()
        self.log_ss(self.ss_surf, self.r_surf)

    def log_ss(self, ss=None, r=None):
        if ss is not None:
            log_format = "{0:20.10E} | {1:20.10E} | {2:20.10E} | {3:20.10E} | {4:20.10E}"
            print(log_format.format(r, *ss))
        else:
            print("{0:20} | {1:20} | {2:20} | {3:20} | {4:20}".format(*readable_strings))

    def log_raw(self, a=0, b=0):
        size = self.ss_profile.shape[1]
        if a>0:
            print("Printing first {0} data entires.".format(a))
            self.log_ss()
            for i in np.arange(0, min(size, a), 1):
                self.log_ss(self.ss_profile[:, i], self.r_profile[i])
        if b>0:
            print("Printing last {0} data entires.".format(b))
            self.log_ss()
            for i in np.arange(max(min(size, size-b), 0), size, 1):
                self.log_ss(self.ss_profile[:, i], self.r_profile[i])
