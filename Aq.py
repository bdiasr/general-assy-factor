import scipy 
import numpy as np
import sympy as sym
import math

from scipy import special
import scipy.integrate as integrate

import scipy.special as s
import matplotlib.pyplot as plt
plt.rcParams['text.usetex'] = True

import os

import torch
from torchquad import Simpson, set_up_backend, Boole, MonteCarlo, VEGAS, Trapezoid
import inspect
import torch.nn as nn

from cte import *


device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
set_up_backend("torch", data_type="float64")
num_gpus = torch.cuda.device_count()
torch.set_printoptions(precision=20)
dimension = 1
simp = Simpson()
integration_domain = [[0, L]] * dimension
Qs = np.arange(-75, 76)
points = 130




def Fz(z):
    z_tensor = torch.as_tensor(z)
    L_tensor = torch.as_tensor(L)

    condition = (z_tensor >= L_tensor / 12) & (z_tensor <= 11 * L_tensor / 12)
    expoente = -5 * ((z_tensor - (0.5 * L_tensor)) ** 2 * 1 / L_tensor ** 2)

    value = torch.where(condition,
                        torch.exp(expoente) * torch.cos(6 * torch.pi * z_tensor / L_tensor),
                        torch.tensor(0))
    return value


def Aq(q, L):

    frac = 1 / L 
    integration_domain = [[0, L]]
    
    integrand = lambda z: (Fz(z) * 
                           torch.exp(torch.as_tensor(1j * z * (2 * torch.pi * q) / L)))
    
    real_integrand = lambda x: integrand(x).real
    complex_integrand = lambda x: integrand(x).imag
    
    real_integral = simp.integrate((real_integrand), dim=1, N=301, integration_domain=integration_domain).item()
    complex_integral = simp.integrate((complex_integrand), dim=1, N=301, integration_domain=integration_domain).item()
    
    if complex_integral == 0.0:
        integral_val = real_integral * frac
    
    else:
        integral_val = (real_integral + 1j*complex_integral) * frac
        
    return integral_val

'''
aq_normal = []
for q in Qs:
    aq_normal.append(np.log(np.abs(Aq(q, L))))

plt.figure(figsize=[7,5])
plt.plot(Qs, aq_normal, 'r.-')
plt.ylim(-15, 0)
plt.xlim(-75,76)
plt.ylabel(r'$A_q$', fontsize=15)
plt.xlabel(r'$q$', fontsize=15)
plt.grid()
plt.show()
'''