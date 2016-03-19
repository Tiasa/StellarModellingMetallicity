from __future__ import division, print_function
from constants import *

def getmu(X, Y):
    """
    Get the mean molecular weight from X, Y
    """
    Z = 1 - X - Y
    isfractional(X, "X")
    isfractional(Y, "Y")
    isfractional(Z, "Z")

    return 1 / (2 * X + 0.75 * Y + 0.5 * Z)

def isfractional(a, A):
    if not (0 < a < 1):
       raise Exception("{A} should be between 0 and 1 (is {a}).".format(A=A, a=a))

print(getmu(0.7, 0.2))
print(getmu(0.7, 0.4))