import numpy as np 
from mpmath import fac, exp, sin, cos, besselj, legenp, sqrt, conj, pi

lamb = 1064*10**(-9)
#lamb = 9.58*10**(-6)
k = (2 * np.pi)/lamb
L = 400*10**(-6)
Q = 0.8*k

M = (1.57 - 0.038j)
M1 = (1.57 - 0.38j)
M2 = 1.57 - 0.01*1j
M3 = 1.57 - 1j

M2_laop = (1.57 - 0.19j)
M3_laop = (1.57 - 0.95j)

mi = 1
#alpha = 10
eta_r = 1/M


def k_rho(k, alpha):
    return (k*np.sin(alpha))

def k_z(k, alpha):
    return (k*np.cos(alpha))

def Z0(k, alpha, x_0, y_0):
    rho_0 = (x_0**2 + y_0**2)**(1/2)
    return np.float64(rho_0 * k_rho(k, alpha))

def epsilonR2(m, ur):
    epsilon = (m**2) / ur
    return - epsilon.imag

def Eta_r(m, ur):
    epsilon = (m**2) / ur
    divUrEpsilon = ur/epsilon 
    return divUrEpsilon ** (1/2)  
