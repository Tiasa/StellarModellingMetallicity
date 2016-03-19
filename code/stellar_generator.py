import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc

# Computer modern fonts
rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)

def test_fnc(a):
    print a

test_fnc("Stars are cool!")

