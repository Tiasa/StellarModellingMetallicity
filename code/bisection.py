from __future__ import division, print_function
import numpy as np
from progress import printProgress

def bisection(f, a, b, tol):
    n_max = np.ceil(np.log2(abs(b-a) / tol))

    n = 1
    printProgress(0, n_max, "Bisection")
    f_a = f(a)
    f_b = f(b)
    best_c = None
    best_f = None
    while n <= n_max:
        c = (a + b)/2
        f_c = f(c)
        if best_c is None or abs(f_c) < abs(best_f):
            best_c = c
            best_f = f_c

        # print("Error {0}".format(f_c))
        # print("rho_c {0}".format(c))
        if (f_c == 0 or (b-a)/2 < tol):
            printProgress(n_max, n_max, "Complete")
            return best_c
        printProgress(n, n_max, "Bisection")
        n = n+1
        if (f_c * f_a > 0):
            a, f_a = c, f_c
        else:
            b, f_b = c, f_c
    raise Exception("n_max was not large enough.")