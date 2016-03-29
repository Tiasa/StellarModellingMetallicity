from __future__ import division, print_function
from multiprocessing import Process, Queue
import multiprocessing
from constants import *
from composition import Composition
import matplotlib.pyplot as plt
from matplotlib import rc
from stellar_generator import Star
from dot_dict import DotDict
from progress import printProgress
from timing_profiler import timing
# Computer modern fonts
rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)

N_CORES = multiprocessing.cpu_count()
PROCESSES_PER_CORE = 2
LOG = True

star_attributes_to_pickle = [ # Must be serializable attributes
    'i_surf',
    'ss_profile',
    'r_profile',
    'ss_surf',
    'r_surf',
    'density_c',
    'temp_c',
    'composition',
    'delta_tau_surf',
    'lumin_surf_bb',
    'lumin_surf_rkf',
    'temp_surf',
    'lumin_surf',
    'mass_surf',
    'density',
    'data_size',
]

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

class MainSequence():

    def __init__(self, min_core_temp, max_core_temp,composition,num_stars):
        self.min_core_temp = min_core_temp
        self.max_core_temp = max_core_temp
        self.num_stars = num_stars
        self.composition = composition
        self.__solved = False

    def star_worker(self, temp_c_vals, out_queue):
        for temp_c in temp_c_vals:
            star = Star(temp_c = temp_c, composition=self.composition)
            star.solve()
            star_pickle = {}
            for attribute in star_attributes_to_pickle:
                star_pickle[attribute] = getattr(star, attribute)
            out_queue.put(star_pickle)
            if LOG: printProgress(out_queue.qsize(), self.num_stars, "Workers")

    @timing
    def solve_stars(self):
        stars = []
        out_queue = Queue()
        num_processes = N_CORES * PROCESSES_PER_CORE
        stars_per_process = np.ceil(self.num_stars / num_processes)
        temp_c_space = np.linspace(start=self.min_core_temp,stop=self.max_core_temp,num=self.num_stars)
        # print("Solving {n} stars with temp_c_space:".format(n=self.num_stars))
        processes = []
        i = 0
        if LOG: printProgress(0, self.num_stars, "Workers")
        while i < self.num_stars:
            remaining_stars = self.num_stars - i
            batch = min(remaining_stars, stars_per_process)
            process = Process(target=self.star_worker, args=(temp_c_space[i:i + batch], out_queue))
            processes.append(process)
            process.start()
            i += batch

        for i in xrange(self.num_stars):
            star = DotDict(out_queue.get())
            stars.append(star)

        for process in processes:
            process.join()

        if LOG: printProgress(self.num_stars, self.num_stars, "Workers")

        self.stars = stars
        self.__solved = True

    def plot(self):
        assert self.__solved, "Stars have not be solved yet."

        padding = 0.2

        temp_surf = np.array([star.temp_surf for star in self.stars])
        lumin_surf = np.array([star.lumin_surf for star in self.stars]) / L_s
        mass_surf = np.array([star.mass_surf for star in self.stars]) / M_s
        r_surf = np.array([star.r_surf for star in self.stars]) / R_s
        mass_space = np.linspace(np.min(mass_surf)*(1 - padding), np.max(mass_surf)*(1 + padding), 300)

        # Main Sequence
        plt.figure()
        plt.title(r"Main Sequence")
        plt.xlabel(r"Temperature (K)")
        plt.ylabel(r"$L/L_{\odot}$")
        plt.plot(temp_surf,lumin_surf, "x")
        plt.gca().invert_xaxis()
        plt.gca().set_yscale("log")
        plt.gca().set_xscale("log")
        plt.savefig("../figures/main_sequence_{0}_stars.pdf".format(self.num_stars), format="pdf")
        plt.show()

        # L/L_sun as a function of M/M_sun
        plt.figure()
        plt.title("$L/L_{\odot}$ as a function of $M/M_{\odot}$")
        plt.xlabel(r"$M/M_{\odot}$")
        plt.ylabel(r"$L/L_{\odot}$")
        plt.plot(mass_space,vlme(mass_space),"r--")
        plt.plot(mass_surf,lumin_surf,"bx")
        plt.gca().set_yscale("log")
        plt.gca().set_xscale("log")
        plt.savefig("../figures/LvM_{0}_stars.pdf".format(self.num_stars), format="pdf")
        plt.show()

        ## R/R_sun as a function of M/M_sun
        plt.figure()
        plt.title("$R/R_{\odot}$ as a function of $M/M_{\odot}$")
        plt.xlabel(r"$M/M_{\odot}$")
        plt.ylabel(r"$R/R_{\odot}$")
        plt.plot(mass_space,vrme(mass_space),"r--")
        plt.plot(mass_surf, r_surf,"bx")
        plt.gca().set_yscale("log")
        plt.gca().set_xscale("log")
        plt.savefig("../figures/RvM_{0}_stars.pdf".format(self.num_stars), format="pdf")
        plt.show()

if __name__ == "__main__":
    LOG_BS = False
    LOG_SG = False
    main_seq = MainSequence(min_core_temp=5e6, max_core_temp=3.5e7, composition=Composition.fromXY(0.73,0.25), num_stars=5)
    main_seq.solve_stars()
    main_seq.plot()
