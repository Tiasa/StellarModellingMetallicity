from __future__ import division, print_function

from matplotlib import rc
import matplotlib.pyplot as plt
from stellar_generator import *

# Computer modern fonts
rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)

def plot_star(star):
    if not star.is_solved:
        star.solve()

    i_surf = star.i_surf
    i_count = len(star.r_profile)
    lumin_surf = star.lumin_surf
    temp_surf = star.temp_surf
    r_surf = star.r_surf
    mass_surf = star.mass_surf
    lumin_p = star.ss_profile[lumin, :]
    mass_p = star.ss_profile[mass, :]
    density_p = star.ss_profile[density, :]
    temp_p = star.ss_profile[temp, :]
    central = star.ss_profile[:, 0]
    r_0 = star.r_profile[0]
    pressure_c = star.pressure(central, r_0)
    pressure_p = np.zeros([3, i_count])
    r = star.r_profile

    n_r = r / r_surf
    n_lumin = lumin_p / lumin_surf
    n_mass = mass_p / mass_surf
    n_temp = temp_p / star.temp_c
    n_density = density_p / star.density_c
    n_ss = [n_density, n_temp, n_lumin, n_mass]

    plt.figure()
    plt.title(r"Normalized Stellar State")
    plt.xlabel(r"Radius ($r/R_*$)")
    plt.ylabel(r"$\rho/\rho_c, L/L_*, T/T_c, M/M_*$")
    n_ss_labels = [r"$\rho$", r"$T$", r"$L$", r"$M$"]
    plots = [plt.plot(n_r, n_ss_i)[0] for n_ss_i in n_ss]
    plt.axis([0,n_r[-1],0, 1.05])
    plt.legend(plots, n_ss_labels, loc="best")
    plt.gca().set_autoscale_on(False)
    plt.axvspan(n_r[i_surf], n_r[-1], color='lightgray', alpha=0.5)
    plt.show()

def plot_step_sizes(star):
    if not star.is_solved:
        star.solve()

    steps = star.r_profile[1:] - star.r_profile[:-1]
    plt.figure()
    plt.title("Step Size Required")
    plt.gca().set_yscale("log")
    plt.plot(range(len(steps)), steps)
    plt.show()


# test_star = Star(temp_c = 1.5e7, density_c=1.6e5, composition=Composition.fromXY(0.69, 0.29))
# test_star = Star(temp_c = 3e7, composition=Composition.fromXY(0.73, 0.25))
# test_star = Star(temp_c = 1.2e10, composition=Composition.fromXY(0.73, 0.25))
# test_star = Star(temp_c = 1e6, composition=Composition.fromXY(0.73, 0.25))
# test_star = Star(temp_c = 3.5e7, composition=Composition.fromXY(0.5, 0.1))
# test_star = Star(temp_c = 8.23e6, composition=Composition.fromXY(0.73, 0.25))
test_star = Star(temp_c = 1e6, composition=Composition.fromXY(0.73, 0.25))

test_star.solve()
# test_star.log_raw(b=20)
test_star.log_solved_properties()

# plot_step_sizes(test_star)
plot_star(test_star)