from __future__ import division, print_function

from matplotlib import rc
import matplotlib.pyplot as plt

# Computer modern fonts
rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)

def plotStar(star):

    def plot(self):
        if not self.__solved:
            self.solve()

        lumin_surf = self.lumin_surf
        temp_surf = self.temp_surf
        r_surf = self.r_surf
        mass_surf = self.mass_surf
        lumin_p = self.ss_profile[lumin, :]
        mass_p = self.ss_profile[mass, :]
        density_p = self.ss_profile[density, :]
        temp_p = self.ss_profile[temp, :]
        r = self.r_profile

        n_r = r / r_surf
        n_lumin = lumin_p / lumin_surf
        n_mass = mass_p / mass_surf
        n_temp = temp_p / temp

        for i in range(self.ss_profile.shape[0]):
            plt.figure()
            plt.title("{0} vs. Radius".format(readable_strings[i+1]))
            plt.xlabel("Radius")
            plt.ylabel(readable_strings[i+1])
            plt.plot(self.r_profile, self.ss_profile[i, :])
            plt.show()

    def plot_step_sizes(self):
        if not self.__solved:
            self.solve()

        steps = self.r_profile[1:] - self.r_profile[:-1]
        plt.figure()
        plt.title("Step Size Required")
        plt.gca().set_yscale("log")
        plt.plot(range(len(steps)), steps)
        plt.show()
