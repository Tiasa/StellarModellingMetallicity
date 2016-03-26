from __future__ import division, print_function

class Composition():

    def __init__(self, X, Y, Z):
        isfractional(X, "X")
        isfractional(Y, "Y")
        isfractional(Z, "Z")
        assert (X + Y + Z == 1)
        self.X = X
        self.Y = Y
        self.Z = Z
        self.mu = 1 / (2 * X + 0.75 * Y + 0.5 * Z)

    @staticmethod
    def fromXY(X,Y):
        return Composition(X, Y, 1 - X - Y)

    @staticmethod
    def fromYZ(Y,Z):
        return Composition(1 - Y - Z, Y, Z)

    @staticmethod
    def fromZX(Z,X):
        return Composition(X, 1 - X - Z, Z)

def isfractional(a, A):
    if not (0 < a < 1):
       raise Exception("{A} should be between 0 and 1 (is {a}).".format(A=A, a=a))

# Tests
# print(Composition.fromXY(0.4,0.3).mu)
# print(Composition.fromXY(0.4,0.9).mu)
