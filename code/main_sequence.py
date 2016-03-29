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
PROCESSES_PER_CORE = 1
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
    'density_surf',
    'data_size',
]

class MainSequence():

    def __init__(self, min_core_temp, max_core_temp,composition,num_stars):
        self.min_core_temp = min_core_temp
        self.max_core_temp = max_core_temp
        self.num_stars = num_stars
        self.composition = composition
        self.solved = False

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
        self.temp_surf = np.array([star.temp_surf for star in self.stars])
        self.lumin_surf = np.array([star.lumin_surf for star in self.stars])
        self.n_lumin_surf = self.lumin_surf / L_s
        self.mass_surf = np.array([star.mass_surf for star in self.stars])
        self.n_mass_surf = self.mass_surf / M_s
        self.r_surf = np.array([star.r_surf for star in self.stars])
        self.n_r_surf = self.r_surf / R_s
        self.solved = True
