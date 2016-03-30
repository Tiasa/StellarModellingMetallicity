from __future__ import division, print_function

from matplotlib import rc
import os
import matplotlib.pyplot as plt
from stellar_generator import *
from constants import gamma,mach_ep
from composition import Composition
from main_sequence import MainSequence
from where_positive import where_positive
import math
import os

# Computer modern fonts
rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)

def plot_star(star):
    if not star.is_solved:
        star.solve()

    star_dir_name = "../figures/star_comp-{comp}_Tc-{temp_c}".format(comp = star.composition.file_string,temp_c=star.temp_c)
    # Previous file name
    #star_file_name = "../figures/{prefix}_star_comp-{comp}_Tc-{temp_c}.pdf".format(prefix = "{prefix}", comp = star.composition.file_string, temp_c = star.temp_c)
    star_file_name = (star_dir_name+"/{prefix}.png").format(prefix="{prefix}")
    if not os.path.exists(os.path.dirname(star_file_name)):
        try:
            os.makedirs(os.path.dirname(star_file_name))
        except OSError as exc: # Guard against race condition
            if exc.errno != errno.EEXIST:
                raise
    i_surf = star.i_surf
    i_count = len(star.r_profile)
    lumin_surf = star.lumin_surf
    temp_surf = star.temp_surf
    r_surf = star.r_surf
    mass_surf = star.mass_surf
    lumin_p = star.ss_profile[lumin, :]
    mass_p = star.ss_profile[mass, :]
    density_p = star.ss_profile[density, :]
    density_c = star.density_c
    temp_p = star.ss_profile[temp, :]
    central = star.ss_profile[:, 0]
    r = star.r_profile
    r_0 = star.r_profile[0]
    pressure_c = star.pressure(central, r_0)
    pressure_p = np.zeros([4, i_count])
    opacity_p = np.zeros([4, i_count])
    dL_dr_p = np.zeros([3, i_count])
    is_convective = np.zeros(i_count)
    radiative = np.zeros(i_count)
    convective = np.zeros(i_count)
    pressure_grad = np.zeros(i_count)
    for i in xrange(i_count):
        ss_i = star.ss_profile[:, i]
        r_i = star.r_profile[i]
        partial_pressure = star.partial_pressure(ss_i, r_i)
        pressure_grad[i] = star.dP_dr(ss_i, r_i)
        pressure_p[0, i] = sum(partial_pressure)
        pressure_p[1:4, i] = [p for p in partial_pressure]

        partial_opacity = star.partial_opacity(ss_i, r_i)
        opacity_p[0, i] = star.opacity(ss_i, r_i)
        opacity_p[1:4, i] = [k for k in partial_opacity]

        partial_dL_dr = star.partial_dL_dr(ss_i, r_i)
        dL_dr_p[0, i] = star.dL_dr(ss_i, r_i)
        dL_dr_p[1:3, i] = [k for k in partial_dL_dr]

        partital_energy_trans = star.partial_dT_dr(ss_i, r_i)
        radiative[i], convective[i] = partital_energy_trans
        if (radiative[i] > convective[i]):
            is_convective[i] = 1

    convective_regions = where_positive(is_convective)

    n_r = r / r_surf
    n_lumin = lumin_p / lumin_surf
    n_mass = mass_p / mass_surf
    n_temp = temp_p / star.temp_c
    n_density = density_p / star.density_c
    n_pressure = pressure_p / pressure_c
    n_opacity = np.log10(opacity_p)
    n_dL_dr = dL_dr_p * r_surf / lumin_surf
    n_ss = [n_density, n_temp, n_lumin, n_mass]

    # dlogP_dlogT = - (temp_p / pressure_p[0]) * pressure_grad / np.minimum(radiative, convective)
    # dlogP_dlogT = - (temp_p / pressure_p[0]) * pressure_grad / np.minimum(radiative, convective)
    logP = np.log(pressure_p[0])[:i_surf]
    logT = np.log(temp_p)[:i_surf]
    dlogP = logP[1:] - logP[:-1]
    dlogT = logT[1:] - logT[:-1]

    dlogP_dlogT = dlogP/(dlogT+mach_ep)

    # Plot for stellar state values
    plt.figure()
    plt.title(r"Normalized Stellar State")
    plt.xlabel(r"Radius ($r/R_*$)")
    plt.ylabel(r"$\rho/\rho_c, L/L_*, T/T_c, M/M_*$")
    n_ss_labels = [r"$\rho$", r"$T$", r"$L$", r"$M$"]
    plots = [plt.plot(n_r, n_ss_i)[0] for n_ss_i in n_ss]
    plt.axis([0,n_r[-1],0, 1.05])
    plt.legend(plots, n_ss_labels, loc="best")
    plt.gca().set_autoscale_on(False)
    for region in convective_regions:
        plt.axvspan(n_r[region[0]], n_r[region[1]], color='gray', alpha=0.4)
    plt.savefig(star_file_name.format(prefix="stellar_state"), format="pdf")
    # plt.show()

    # Plotting pressure decomposition
    plt.figure()
    plt.title(r"Partial Pressures")
    plt.xlabel(r"Radius ($r/R_*$)")
    plt.ylabel(r"$P/P_c$")
    n_pp_labels = [r"$P_{\mathrm{tot}}$", r"$P_{\mathrm{deg}}$", r"$P_{\mathrm{gas}}$", r"P_{\gamma}"]
    plots = [plt.plot(n_r, n_pressure_i)[0] for n_pressure_i in n_pressure]
    plt.axis([0,n_r[-1],0, 1.05])
    plt.legend(plots, n_pp_labels, loc="best")
    plt.gca().set_autoscale_on(False)
    for region in convective_regions:
        plt.axvspan(n_r[region[0]], n_r[region[1]], color='gray', alpha=0.4)
    plt.savefig(star_file_name.format(prefix="partial_pressure"), format="pdf")
    # plt.show()

    # Plotting luminosity decomposition
    plt.figure()
    plt.title(r"Partial Luminosities")
    plt.xlabel(r"Radius ($r/R_*$)")
    plt.ylabel(r"$\mathrm{d}L/\mathrm{d}r$ ($L_*/R_*$)")
    n_pl_labels = [r"$\mathrm{d}L/\mathrm{d}r$", r"$\mathrm{d}L_{\mathrm{PP}}/\mathrm{d}r$", r"$\mathrm{d}L_{\mathrm{CNO}}/\mathrm{d}r$"]
    plots = [plt.plot(n_r, n_dL_dr_i)[0] for n_dL_dr_i in n_dL_dr]
    plots[1].set_dashes([4,4])
    plots[2].set_dashes([8,4,2,4])
    plt.axis([0,n_r[-1],0,max(n_dL_dr[0])*1.5])
    plt.legend(plots, n_pl_labels, loc="best")
    plt.gca().set_autoscale_on(False)
    for region in convective_regions:
        plt.axvspan(n_r[region[0]], n_r[region[1]], color='gray', alpha=0.4)
    plt.savefig(star_file_name.format(prefix="partial_lumin"), format="pdf")
    # plt.show()

    # Plotting opacity decomposition
    plt.figure()
    plt.title(r"Partial Opacities")
    plt.xlabel(r"Radius ($r/R_*$)")
    plt.ylabel(r"$\log_{10}(\kappa)$")
    n_pk_labels = [r"$\kappa_{\mathrm{tot}}$", r"$\kappa_{\mathrm{es}}$", r"$\kappa_{\mathrm{ff}}$", r"$\kappa_{\mathrm{H}^-}$"]
    plots = [plt.plot(n_r, n_opacity_i)[0] for n_opacity_i in n_opacity]
    plt.axis([0,n_r[-1], np.min(n_opacity), np.max(n_opacity)])
    plt.legend(plots, n_pk_labels, loc="best")
    plt.gca().set_autoscale_on(False)
    for region in convective_regions:
        plt.axvspan(n_r[region[0]], n_r[region[1]], color='gray', alpha=0.4)
    plt.savefig(star_file_name.format(prefix="partial_opacity"), format="pdf")
    # plt.show()

    # Plotting logPlogT
    plt.figure()
    plt.title(r"Convective Stability")
    plt.xlabel(r"Radius ($r/R_*$)")
    plt.ylabel(r"$\mathrm{d}\log P/\mathrm{d}\log T$")
    n_lPlT_labels = [r"Convective boundary $(1 - 1/\gamma)^{-1}$", r"Star"]
    plots = []
    log_points = len(dlogP_dlogT)
    boundary = np.zeros(log_points)
    boundary.fill(1/(1-1/gamma))
    plots.append(plt.plot(n_r[:log_points], boundary)[0])
    plots.append(plt.plot(n_r[:log_points], dlogP_dlogT)[0])
    plt.axis([0,n_r[log_points], 0, np.max(dlogP_dlogT) * 1.3])
    plt.legend(plots, n_lPlT_labels, loc="best")
    plt.gca().set_autoscale_on(False)
    for region in convective_regions:
        plt.axvspan(n_r[region[0]], n_r[region[1]], color='gray', alpha=0.4)
    plt.savefig(star_file_name.format(prefix="dlogP_dlogT"), format="pdf")
    # plt.show()
    ## Saving the star specifics
    ## We are targetting :
    ## 1. Surface Temp
    ## 2. Central density
    ## 3. Central Temperature
    ## 4. Radius
    ## 5. Mass
    ## 6. Luminosity
    f = open(star_dir_name+'/profile.txt', 'w')
    f.write('Surface Temperature = '+ repr(temp_surf) + '\n')
    f.write('Central Density = '+repr(density_c)+'\n')
    f.write('Radius = '+repr(r_surf)+'\n')
    f.write('Mass = '+repr(mass_surf/M_s)+'\n')
    f.write('Luminosity = '+ repr(lumin_surf) +'\n')
    f.close()


def plot_step_sizes(star):
    if not star.is_solved:
        star.solve()

    steps = star.r_profile[1:] - star.r_profile[:-1]
    plt.figure()
    plt.title("Step Size Required")
    plt.gca().set_yscale("log")
    plt.plot(range(len(steps)), steps)
    plt.show()

def lumin_mass_exact(m):
    if m < 0.7:
        return 0.35*m**2.62
    else:
        return 1.02*m**3.92
vlme = np.vectorize(lumin_mass_exact)

def radius_mass_exact(m):
    if m < 1.66:
        return 1.06*m**0.945
    else:
        return 1.33*m**0.555
vrme = np.vectorize(radius_mass_exact)

def plot_main_sequence(v_main_seq):
    for main_seq in v_main_seq:
        if not main_seq.solved:
            main_seq.solve()


    labels = [r"$Z = {Z}$".format(Z=main_seq.composition.Z) for main_seq in v_main_seq]
    num_stars = sum([main_seq.num_stars for main_seq in v_main_seq])
    num_seqs = len(v_main_seq)

    main_seq_folder = "../figures/main_sequence_{0}_stars_{1}_seq/".format(num_stars, num_seqs)
    if not os.path.exists(main_seq_folder):
        os.makedirs(main_seq_folder)

    padding = 0.2
    # Main Sequence
    plt.figure()
    plt.title(r"Main Sequence")
    plt.xlabel(r"Temperature (K)")
    plt.ylabel(r"$L/L_{\odot}$")
    #eps = 1e-20
    #blah = log()
    blah = [math.log(x) for x in main_seq.temp_surf]
    plots = [plt.plot(blah, main_seq.n_lumin_surf, "+")[0] for main_seq in v_main_seq]
    plt.legend(plots, labels, loc="best")
    plt.gca().invert_xaxis()
    plt.gca().set_yscale("log")
    #plt.gca().set_xscale("log")
    plt.savefig(main_seq_folder + "ms.pdf", format="pdf")
    plt.show()

    merge_mass_surf = [main_seq.n_mass_surf for main_seq in v_main_seq]
    mass_space = np.linspace(np.min(merge_mass_surf)*(1 - padding), np.max(merge_mass_surf)*(1 + padding), 300)

    # L/L_sun as a function of M/M_sun
    plt.figure()
    plt.title("$L/L_{\odot}$ as a function of $M/M_{\odot}$")
    plt.xlabel(r"$M/M_{\odot}$")
    plt.ylabel(r"$L/L_{\odot}$")
    plots = [plt.plot(main_seq.n_mass_surf,main_seq.n_lumin_surf,"+")[0] for main_seq in v_main_seq]
    plots.append(plt.plot(mass_space,vlme(mass_space),"r--")[0])
    plt.legend(plots, labels + ["Empirical"], loc="best")
    plt.gca().set_yscale("log")
    plt.gca().set_xscale("log")
    plt.savefig(main_seq_folder + "LvM.pdf", format="pdf")
    plt.show()

    ## R/R_sun as a function of M/M_sun
    plt.figure()
    plt.title("$R/R_{\odot}$ as a function of $M/M_{\odot}$")
    plt.xlabel(r"$M/M_{\odot}$")
    plt.ylabel(r"$R/R_{\odot}$")
    plots = [plt.plot(main_seq.n_mass_surf,main_seq.n_r_surf,"+")[0] for main_seq in v_main_seq]
    plots.append(plt.plot(mass_space,vrme(mass_space),"r--")[0])
    plt.legend(plots, labels + ["Empirical"], loc="best")
    plt.gca().set_yscale("log")
    plt.gca().set_xscale("log")
    plt.savefig(main_seq_folder + "RvM.pdf", format="pdf")
    plt.show()

if __name__ == "__main__":
    # Remember to turn off logging in adaptive_bisection.py
#    composition = [Composition.fromZX(Z, 0.73) for Z in [0.00, 0.01, 0.015, 0.02, 0.03]]
#    v_main_seq = [MainSequence(min_core_temp=5e6, max_core_temp=3.5e7, composition=comp, num_stars=100) for comp in composition]
#    for main_seq in v_main_seq:
#        main_seq.solve_stars()

#    plot_main_sequence(v_main_seq)

    # test_star = Star(temp_c = 1.5e7, density_c=1.6e5, composition=Composition.fromXY(0.69, 0.29))
    # test_star = Star(temp_c = 3e7, composition=Composition.fromXY(0.73, 0.25))
    # test_star = Star(temp_c = 1.2e10, composition=Composition.fromXY(0.73, 0.25))
    # test_star = Star(temp_c = 1e6, composition=Composition.fromXY(0.73, 0.25))
    # test_star = Star(temp_c = 3.5e7, composition=Composition.fromXY(0.5, 0.1))
    # test_star = Star(temp_c = 1e8, composition=Composition.fromXY(0.73, 0.25))
    # test_star = Star(temp_c = 3.5e7, composition=Composition.fromXY(0.73, 0.25))
     test_star = Star(temp_c = 8.23e6, composition=Composition.fromXY(0.65, 0.25))

    # test_star.solve()
    # # # test_star.log_raw(b=20)
    # test_star.log_solved_properties()

   # # # plot_step_sizes(test_star)
     plot_star(test_star)
