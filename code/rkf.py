from __future__ import division, print_function
import numpy as np
import matplotlib.pyplot as plt

# Coefficients used to compute the independent variable argument of f

a2  =   2.500000000000000e-01  #  1/4
a3  =   3.750000000000000e-01  #  3/8
a4  =   9.230769230769231e-01  #  12/13
a5  =   1.000000000000000e+00  #  1
a6  =   5.000000000000000e-01  #  1/2

# Coefficients used to compute the dependent variable argument of f

b21 =   2.500000000000000e-01  #  1/4
b31 =   9.375000000000000e-02  #  3/32
b32 =   2.812500000000000e-01  #  9/32
b41 =   8.793809740555303e-01  #  1932/2197
b42 =  -3.277196176604461e+00  # -7200/2197
b43 =   3.320892125625853e+00  #  7296/2197
b51 =   2.032407407407407e+00  #  439/216
b52 =  -8.000000000000000e+00  # -8
b53 =   7.173489278752436e+00  #  3680/513
b54 =  -2.058966861598441e-01  # -845/4104
b61 =  -2.962962962962963e-01  # -8/27
b62 =   2.000000000000000e+00  #  2
b63 =  -1.381676413255361e+00  # -3544/2565
b64 =   4.529727095516569e-01  #  1859/4104
b65 =  -2.750000000000000e-01  # -11/40

# Coefficients used to compute local truncation error estimate.  These
# come from subtracting a 4th order RK estimate from a 5th order RK
# estimate.

r1  =   2.777777777777778e-03  #  1/360
r3  =  -2.994152046783626e-02  # -128/4275
r4  =  -2.919989367357789e-02  # -2197/75240
r5  =   2.000000000000000e-02  #  1/50
r6  =   3.636363636363636e-02  #  2/55

# Coefficients used to compute 4th order RK estimate

c1  =   1.157407407407407e-01  #  25/216
c3  =   5.489278752436647e-01  #  1408/2565
c4  =   5.353313840155945e-01  #  2197/4104
c5  =  -2.000000000000000e-01  # -1/5

BUFFER = 2**14

def rkf( f, a, x0, tol, stop ):
    """Vectorized Runge-Kutta-Fehlberg method to solve x' = f(x,t) with x(t[0]) = x0.

    USAGE:
        T, X = rkf(f, a, b, x0, tol)

    INPUT:
        f     - an array that is the system of equations dvx/dt = vf(vx,t)
        a     - left-hand endpoint of interval (initial condition is here)
        x0    - initial x value: x0 = x(a)
        tol   - maximum value of local truncation error estimate
        stop  - stopping condition func(i, vx, t)

    OUTPUT:
        T     - NumPy array of independent variable values
        X     - NumPy array of corresponding solution function values

    NOTES:
        This function implements 4th-5th order Runge-Kutta-Fehlberg Method
        to solve the initial value problem

           dx
           -- = f(x,t),     x(a) = x0
           dt

        on the interval [a,...).

        Based on pseudocode presented in "Numerical Analysis", 6th Edition,
        by Burden and Faires, Brooks-Cole, 1997.
    """

    s = len(x0)

    # Set t and x according to initial condition and assume that h begins with resonable value

    t = a
    x = x0
    # h = hmin
    h = 1e-6

    # max_samples = np.ceil(abs((b-a)/hmin))

    # print("Sampling with initial buffer {0}".format(BUFFER))

    # Initialize arrays that will be returned

    T = np.empty( BUFFER )
    X = np.empty( [s, BUFFER] )

    T[0] = t
    # print(x.shape)
    # print(X.shape)
    X[:,0] = x

    i = 0

    while not stop(i,x,t):
        if i + 1 >= len(T):
            # print("Increasing buffer by {0}".format(BUFFER))
            T = np.hstack((T, np.empty(BUFFER)))
            X = np.hstack((X, np.empty( [s, BUFFER] )))

        # Compute values needed to compute truncation error estimate and
        # the 4th order RK estimate.

        try:
            k1 = h * f( x, t )
            k2 = h * f( x + b21 * k1, t + a2 * h )
            k3 = h * f( x + b31 * k1 + b32 * k2, t + a3 * h )
            k4 = h * f( x + b41 * k1 + b42 * k2 + b43 * k3, t + a4 * h )
            k5 = h * f( x + b51 * k1 + b52 * k2 + b53 * k3 + b54 * k4, t + a5 * h )
            k6 = h * f( x + b61 * k1 + b62 * k2 + b63 * k3 + b64 * k4 + b65 * k5, t + a6 * h )
        except ValueError:
            # The step size must be too large and it is interpolating to negative values
            h = h * 0.1
            continue

        # Compute the estimate of the local truncation error.  If it's small
        # enough then we accept this step and save the 4th order estimate.

        r = abs( r1 * k1 + r3 * k3 + r4 * k4 + r5 * k5 + r6 * k6 ) / h
        if len( np.shape( r ) ) > 0:
            r = max( r )
        if r <= tol:
            t = t + h
            x = x + c1 * k1 + c3 * k3 + c4 * k4 + c5 * k5
            T[i] = t
            X[:,i] = x
            i += 1

        # Now compute next step size, and make sure that it is not too big or
        # too small.

        # h = h * min( max( 0.84 * ( tol / (r + np.finfo(float).eps) )**0.25, 0.1 ), 4.0 )
        # h = h * min( 0.84 * ( tol / (r + np.finfo(float).eps) )**0.25, 4.0 )
        h = h * 0.84 * ( tol / (r + np.finfo(float).eps) )**0.25

    # endwhile

    T = T[0:i]
    X = X[:,0:i]

    # print("Truncating buffer to {0} filled samples.".format(i-1))
    return ( T, X )


# Tests
def test_rkf():
    f0 = lambda x, t: x[1]
    f1 = lambda x, t: -x[0]
    f2 = lambda x, t: x[3]
    f3 = lambda x, t: -x[2]+x[1]/(x[4]**2+1)
    f4 = lambda x, t: x[1]

    stop = lambda i, x, t: t > 10*np.pi
    f = lambda x, t: np.array([f(x,t) for f in [f0, f1, f2, f3, f4]])

    T, X = rkf(f, 0, [1, 0, 0, 0, 0], 1e-6, stop)

    plt.figure()
    for i in range(X.shape[0]):
        plt.plot(T, X[i,:])
    plt.show()

# test_rkf()