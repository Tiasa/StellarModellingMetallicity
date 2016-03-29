from __future__ import division, print_function
from multiprocessing import Process, Lock, Array
from constants import *
from composition import Composition
import matplotlib.pyplot as plt
from matplotlib import rc
from stellar_generator import Star
# Computer modern fonts
rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)

# Main sequence class
# USAGE :
#       ms = MainSequence(min_core_temp = $some temp$, max_core_temp = $some temp$, X = $some comp$, Y = $some comp$, num_of_stars = $some num$)
#
# INPUT:
#       You know
#
# OUTPUT:
#       Plot !!!


class MainSequence():
    def __init__(self, min_core_temp, max_core_temp,X,Y,num_of_stars):
        self.min_core_temp = min_core_temp
        self.max_core_temp = max_core_temp
        self.num_of_stars = num_of_stars
        self.composition = Composition.fromXY(X,Y)
        self.lumin = [0] * num_of_stars
        self.lock = Lock()
    def make_entry(self,t,l,m,r,ind,temp,lumin,mass,radius):
        self.lock.acquire()
        temp[ind] = t
        lumin[ind] = l/L_s
        mass[ind] = m/M_s
        radius[ind] = r/R_s
        print("Star {0} : Temp {1} : Lumin {2} : Mass {3} : Radius {4}\n".format(ind, t, lumin[ind], mass[ind] , radius[ind]))
        self.lock.release()
    def make_star(self,t,index,temp,lumin,mass,radius):
        star =  Star(temp_c = t, composition=self.composition)
        star.solve()
        #star.log_solved_properties()
        self.make_entry(star.temp_surf,star.lumin_surf,star.mass_surf,star.r_surf,index,temp,lumin,mass,radius)
    def calculate(self,temp,lumin,mass,radius):
        core_temp = np.linspace(start=self.min_core_temp,stop=self.max_core_temp,num=self.num_of_stars)
        #print(core_temp)
        procs = []
        for idx, t in enumerate(core_temp):
             p = Process(target=self.make_star,args=(t,idx,temp,lumin,mass,radius))
             procs.append(p)
             p.start()
        for p in procs:
            p.join()
    def calcLuminExact(self,m):
        if m < 0.7:
            return 0.35*m**2.62
        else:
            return 1.02*m**3.92
    def calcRadExact(self,m):
        if m < 1.66:
            return 1.06*m**0.945
        else:
            return 1.33*m**0.555
    def plot(self):
        ## For synchromization create process safe array
        temp = Array("d",[0]*self.num_of_stars)
        lumin = Array("d",[0]*self.num_of_stars)
        mass = Array("d",[0]*self.num_of_stars)
        radius = Array("d",[0]*self.num_of_stars)
        self.calculate(temp,lumin,mass,radius)
        ## HR diagram
        tempArr = temp[:] ## Temp t
        luminArr = lumin[:] ## L/L_sun
        plt.figure()
        plt.title(r"Main Sequence")
        plt.xlabel(r"Temparature")
        plt.ylabel(r"$L/L_{\odot}$")
        plt.plot(tempArr,luminArr,"ro")
        plt.gca().invert_xaxis()
        plt.gca().set_yscale("log")
        #plt.gca().set_xscale("log")
        plt.savefig("../figures/main_sequence_{0}_stars.pdf".format(self.num_of_stars), format="pdf")
        plt.show()
        ## L/L as a function of M/M
        massArr = mass[:] # M/M_sun
        luminExact = [self.calcLuminExact(m) for m in massArr]
        plt.figure()
        plt.title("$L/L_{\odot}$ as a function of $M/M_{\odot}$")
        plt.xlabel(r"$M/M_{\odot}$")
        plt.ylabel(r"$L/L_{\odot}$")
        plt.plot(massArr,luminArr,"ro")
        plt.plot(massArr,luminExact,"bo")
        plt.gca().set_yscale("log")
        plt.gca().set_xscale("log")
        plt.savefig("../figures/LAsAFuncOfM_{0}_stars.pdf".format(self.num_of_stars), format="pdf")
        plt.show()
        ## R/R as a function of M/M
        radiusArr = radius[:] ## R/R_sun
        radiusExact = [self.calcRadExact(m) for m in massArr]
        plt.figure()
        plt.title("$R/R_{\odot}$ as a function of $M/M_{\odot}$")
        plt.xlabel(r"$M/M_{\odot}$")
        plt.ylabel(r"$R/R_{\odot}$")
        plt.plot(massArr,radiusArr,"ro")
        plt.plot(massArr,radiusExact,"bo")
        plt.gca().set_yscale("log")
        plt.gca().set_xscale("log")
        plt.savefig("../figures/RAsAFuncOfM_{0}_stars.pdf".format(self.num_of_stars), format="pdf")
        plt.show()
if __name__ == "__main__":
    ms = MainSequence(min_core_temp=5e6,max_core_temp=3.5e7,X=0.73,Y=0.25,num_of_stars=50)
    ms.plot()
