from constants import *
def dPdrho = (T, rho, mu)
	''' partial derivative P rho '''
	return (1/3)*(3*pi**2)**(2/3)*hbar**2*(1/m_e)*(rho/m_p)*(2/3) + k*T/(mu*m_p)
def dPdT = (rho, mu, T)
	''' partial derivative P T '''
	return rho*k/(mu*m_p) + 4/3*a*T**3
