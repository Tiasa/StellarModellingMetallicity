from __future__ import division, print_function
from multiprocessing import Process, Lock
from constants import *
from composition import Composition
from rkf import rkf
from bisection import bisection
import matplotlib.pyplot as plt
from stellar_generator import Star

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
        self.temp = [0] * num_of_stars
        self.lumin = [0] * num_of_stars
    def make_entry(self,lk,temp,lumin,ind):
        lk.acquire()
        self.temp[ind] = temp
        self.lumin[ind] = lumin
        print("Star {0} : Temp {1} : Lumin {2}".format(ind, temp, lumin))
        lk.release()
    def make_star(self,temp,lk,index):
        star =  Star(temp_c = temp, composition=self.composition)
        star.solve()
         self.make_entry(lk,star.temp_surf,star.lumin_surf,index)
    def calculate(self):
        core_temp = np.linspace(start=self.min_core_temp,stop=self.max_core_temp,num=self.num_of_stars)
        ## For synchromization creating a lock
        lock = Lock()
        for idx, t in enumerate(core_temp):
             Process(target=self.make_star,args=(t,lock,idx)).start()
        
    #def plot(self):
        # DO NOTHING FOR NOW
        
ms = MainSequence(min_core_temp=3.5e7,max_core_temp=9e5,X=0.73,Y=0.25,num_of_stars=10)
ms.calculate()
