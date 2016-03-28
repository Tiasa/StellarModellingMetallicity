from __future__ import division, print_function
import numpy as np
from progress import printProgress

eval_tol_max = 0.6
eval_tol_min = 0.01

def tween(i, a, b):
    assert (0 <= i <= 1), "i needs to be normalized"

    return a + i**4 * (b - a)

def adaptive_bisection(f, a, b, precision):
    n_max = np.ceil(np.log2(abs(b-a) / precision))

    n = 1
    printProgress(0, n_max, "Bisection")
    f_a = f(a, eval_tol_max)
    f_b = f(b, eval_tol_max)
    assert f_a * f_b < 0, "No root in range"
    best_c = None
    best_f = None
    while n <= n_max:
        c = (a + b)/2
        f_c = f(c, tween(n/n_max, eval_tol_max, eval_tol_min))

        if best_c is None or (abs(f_c) < abs(best_f)):
            best_c = c
            best_f = f_c

        if (f_c == 0 or (b-a)/2 < precision):
            printProgress(n_max, n_max, "Complete")
            return best_c
        printProgress(n, n_max, "Bisection")
        n = n+1
        if (f_c * f_a > 0):
            a, f_a = c, f_c
        else:
            b, f_b = c, f_c
    raise Exception("n_max was not large enough.")