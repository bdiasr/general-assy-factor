import os
import time

import numpy as np
import sympy as sym

import scipy as sp 
import scipy.integrate as integrate
import scipy.special as s

import cmath
import math

from mpmath import fac, exp, sin, cos, besselj

import matplotlib.pyplot as plt
plt.rcParams['text.usetex'] = True

from legendre import *
from specials import * 
from cte import * 
from Aq import * 


def g_mnTE(m, n, v, alpha1, lamb, rho_0, phi_0, x_0, y_0, z_0):
    
    cosAxicon =  cos(alpha1)
    sinAxicon =  sin(alpha1)
    
    alpha = alpha1

    a1 = 1j*(-1/2)*(1j**(m+1))*(-1)**((m - abs(m)) * 1/2)   
    b = (fac(n - m))/(fac(n + abs(m)))
    c = exp(1j*k*cosAxicon*z_0)
    d = besselj(m-v-1, k*sinAxicon*rho_0)*exp(-1j*(m-v-1)*phi_0)
    e = (tau_mn(m,n, cosAxicon)/cosAxicon) + (m*pi_mn(m, n, alpha))
    f = besselj(m-v+1, k*sinAxicon*rho_0)*exp(-1j*(m-v+1)*phi_0)
    g = (-tau_mn(m, n, cosAxicon)/cosAxicon) + (m*pi_mn(m, n, alpha))
    
    result = ((a1*b*c)*((d*e)-(f * g)))
    
    return result

def g_mnTM(m, n, v, alpha1, lamb, rho_0, phi_0, x_0, y_0, z_0):
    
    cosAxicon =  cos(alpha1)
    sinAxicon =  sin(alpha1)
    
    alpha = alpha1

    a1 = (1/2) * (1j**(m+1)) * (-1)**((m - abs(m)) * (1/2))
    b = (fac(n - m))/(fac(n + abs(m))) 
    c = exp(1j*k*cosAxicon*z_0) 
    d = (besselj(m - v - 1,  k*sinAxicon*rho_0))*exp(-1j*(m - v - 1)*phi_0)
    e = (tau_mn(m, n, cosAxicon)/cosAxicon) + (m*pi_mn(m, n, alpha))
    f = (besselj(m - v + 1,  k*sinAxicon*rho_0))*exp(-1j*(m - v + 1)*phi_0) 
    g = (-tau_mn(m, n, cosAxicon)/cosAxicon) + (m*pi_mn(m, n, alpha)) 
    
    result =  (a1*b*c)*((d*e) + f*g)
    
    return result

'''
axicon = np.pi/360
spoti = ((405/100)/(k * math.sin(axicon)))
'''


'''
print('Versão Mpmath')
start_time = time.time()
print(g_mnTE(-1, 1, 0,  np.pi/360, lamb, 0.1*spoti,0, 0, 0, 0))
end_time = time.time()
print('Tempo de execução (Mpmath): {} segundos'.format(end_time - start_time))

print('Versão Mpmath')
start_time = time.time()
print(g_mnTM(-1, 1, 0,  np.pi/360, lamb, 0.1*spoti,0, 0, 0, 0))
end_time = time.time()
print('Tempo de execução (Mpmath): {} segundos'.format(end_time - start_time))
'''