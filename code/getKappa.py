import constants
"""
USAGE:
      k = kappa(X,Y,Z,rho,T)
INPUT:
      X , Y , Z : Coefficients
      rho : density
      T : Temparature
OUTPUT:
      kappa : Optical density depending on density and temparature
"""
def kappa(X,Y,Z,rho,T):
    
    ## Electron Scattering Opacity
    kes = 0.2 * (1+X)
    ## Free-free Opacity
    kff = 1e24 * (Z+0.0001) * (rho**0.7)*(T**-3.5)
    ## Hydrogen opacity
    kH = 2.5e-32 * (Z/0.02) * (rho**0.5) * (T**9)
    ## Kappa
    k = (1/kH)+(1/max(kff,kes))
    return k
