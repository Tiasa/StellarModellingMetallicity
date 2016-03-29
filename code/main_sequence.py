from __future__ import division, print_function
from multiprocessing import Process, Queue
from constants import *
from composition import Composition
from rkf import rkf
from bisection import bisection
import matplotlib.pyplot as plt
from matplotlib import rc
from stellar_generator import Star
from dot_dict import DotDict

# Computer modern fonts
rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)

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

class MainSequence():

    def __init__(self, min_core_temp, max_core_temp,composition,num_stars):
        self.min_core_temp = min_core_temp
        self.max_core_temp = max_core_temp
        self.num_stars = num_stars
        self.composition = composition
        self.__solved = False

    def star_worker(self, temp_c, out_queue):
        star = Star(temp_c = temp_c, composition=self.composition)
        star.solve()
        star_pickle = {}
        for attribute in star_attributes_to_pickle:
            star_pickle[attribute] = getattr(star, attribute)
        out_queue.put(star_pickle)

    def solve_stars(self):
        stars = []
        out_queue = Queue()
        temp_c_space = np.linspace(start=self.min_core_temp,stop=self.max_core_temp,num=self.num_stars)
        print("Solving {n} stars with temp_c_space:".format(n=self.num_stars))
        print(temp_c_space)
        processes = []
        for temp_c in temp_c_space:
            process = Process(target=self.star_worker, args=(temp_c, out_queue))
            processes.append(process)
            process.start()

        for i in xrange(self.num_stars):
            star = DotDict(out_queue.get())
            stars.append(star)

        for process in processes:
            process.join()

        self.stars = stars
        self.__solved = True

    def plot(self):
        assert self.__solved, "Stars have not be solved yet."

        temp_surf = np.array([star.temp_surf for star in self.stars])
        lumin_surf = np.array([star.lumin_surf for star in self.stars]) / L_s

        plt.figure()
        plt.title(r"Main Sequence")
        plt.xlabel(r"Temparature")
        plt.ylabel(r"$L/L_{\odot}$")
        plt.plot(temp_surf,lumin_surf, "x")
        plt.gca().invert_xaxis()
        plt.gca().set_yscale("log")
        plt.gca().set_xscale("log")
        plt.savefig("../figures/main_sequence_{0}_stars.pdf".format(self.num_stars), format="pdf")
        plt.show()

if __name__ == "__main__":
    main_seq = MainSequence(min_core_temp=5e6, max_core_temp=3.5e7, composition=Composition.fromXY(0.73,0.25), num_stars=10)
    main_seq.solve_stars()
    main_seq.plot()
    # ms.plot()
