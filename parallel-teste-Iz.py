import os
import time

import numpy as np

import math

from mpmath import fac, exp, sin, cos, besselj, legenp, sqrt, conj, pi

import matplotlib.pyplot as plt
plt.rcParams['text.usetex'] = True

from legendre import *
from specials import * 
from cte import * 
from Aq import * 
from beamShape import *
from coefficients import * 



ceillingX = lambda x : math.ceil(x + (4.05 * (x ** (1/3))) + 2)

def I_z(M, x, epislon, v, alpha, rho_0, z_0):
    
    a = x/k
    
    eta_r = Eta_r(M, 1)

    phi_0 = 0
    
    const = (6* epislon) / ((abs(M)**2) * x**3)
    
    soma = []
    total = 0
    nMax = ceillingX(x)
    x_0 = 0
    y_0 = 0
    
    for n in range(1, nMax+1):
        for m in range(-n, n+1):

            a1 = D_mn(m, n, M, x, v, alpha, rho_0, phi_0, x_0, y_0, z_0)
            b = np.conjugate(S_n(M, n, x)) + (( (n + 1)/M) * r_n(M, n+1, x))

            c = E_mn(m, n, M, x, v, alpha, rho_0, phi_0, x_0, y_0, z_0)

            d = (-S_n(M, n, x) + (((n + 1)/M)*r_n(M, n, x)) )

            e = (1j*F_mn(m, n, M, x, v, alpha, rho_0, phi_0, x_0, y_0, z_0)) 
            sn = S_n(M, n, x)

            soma = (a1*(b)) + (c*(d)) + (e*sn)

            total += soma
            
    return const * total.imag

def iz_serial():
    inicio = time.perf_counter()

    I_z(M, 1, epsilonR2(M, 1), 0, math.radians(5), 0, 0)

    fim = time.perf_counter()

    total = fim - inicio
    print(f"tempo total serial: {total} segundos")


if __name__ == '__main__':

    I_z(M, 1, epsilonR2(M, 1), 0, math.radians(5), 0, 0)