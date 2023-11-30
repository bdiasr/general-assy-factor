
import numpy as np
import scipy.integrate as integrate
import scipy.special as s
import time
import numpy as np
import scipy.integrate as integrate
import scipy.special as s
import math
import torch
from torchquad import Simpson, set_up_backend

from cte import *

import warnings
warnings.filterwarnings(action='ignore', category=UserWarning)

device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
set_up_backend("torch", data_type="float64")
torch.set_printoptions(precision=20)
dimension = 1
simp = Simpson()
integration_domain = [[0, L]] * dimension


def F_z1_pytorch(z, L):
    z_tensor = torch.as_tensor(z)
    L_tensor = torch.as_tensor(L)

    condition = (z_tensor >= L_tensor / 12) & (z_tensor <= 11 * L_tensor / 12)
    expoente = -5 * ((z_tensor - (0.5 * L_tensor)) ** 2 * 1 / L_tensor ** 2)

    value = torch.where(condition,
                        torch.exp(expoente) * torch.cos(6 * torch.pi * z_tensor / L_tensor),
                        torch.tensor(0))
    return(value)


def complex_integrate_pytorch(L, q, a, b):

    integration_domain = [[a, b]]
    
    integrand = lambda z: (F_z1_pytorch(z, L) * 
                           torch.exp(torch.as_tensor(1j*z*(2*torch.pi*q)/L)))
    
    real_integrand = lambda x: integrand(x).real
    complex_integrand = lambda x: integrand(x).imag
    
    real_integral = simp.integrate((real_integrand), dim=1, N=201, integration_domain=integration_domain).item()
    complex_integral = simp.integrate((complex_integrand), dim=1, N=201, integration_domain=integration_domain).item()
    
    if complex_integral == 0.0:
        integral_val = real_integral
    else:
        integral_val = (real_integral + 1j*complex_integral)
        
    return integral_val


def Aq_pytorch(q, L):

    frac = 1/L
    integ = complex_integrate_pytorch(L, q, 0, L)

    return (frac*integ)

def Psi_pytorch(rho, z):
    Q = torch.tensor(0.8*k)
    L = torch.tensor(400*10**(-6))
    soma = []
    total = 0
    qs = np.arange(-75, 76)
    
    for q in qs:

        k_zq = Q + 2 * torch.pi * q/L
        k_pq = torch.sqrt(torch.tensor(k**2 - k_zq**2))
        a = Aq_pytorch(q, L)
        j0 = (torch.special.bessel_j0((k_pq * rho)))
        exponencial = torch.exp((-1j * k_zq * z))
        soma = a * j0 * exponencial
      
        total += soma

    return total


def F_z(z, L):

    if (z >= L/12 and z <= 11*L/12):
        expoente = -5*((z - 0.5*L)**2*1/L**2)
        value = np.exp(expoente) * np.cos(6*np.pi*z/L)
    else:
        value = 0

    return(value)

def complex_integrate(L, q, a, b):

    func = lambda iz: F_z(iz, L) * np.exp( 1j * ((2*(np.pi)*q)/L) * iz)

    real_func = lambda x: func(x).real
    imag_func = lambda x: func(x).imag

    real_integral = integrate.quad(real_func, a, b)
    imag_integral = integrate.quad(imag_func, a, b)

    if imag_integral[0] == 0.0:
        return real_integral[0]

    return real_integral[0] + 1j*imag_integral[0]

def Aq(q, L):

    frac = 1/L
    integ = complex_integrate(L, q, 0, L)

    return frac*integ

def psi(rho, z):

    Q = 0.8*k
    L = 400*10**(-6)
    soma = []
    total = 0
    N = 75
    qs = np.arange(-75, 75)

    for q in qs:

        k_zq = Q + 2*np.pi*q/L
        k_pq = math.sqrt(k**2 - k_zq**2)
        a = Aq(q, L)
        j0 = s.j0(k_pq *rho)
        exponencial = np.exp(-1j * k_zq * z)

        soma = a * j0 * exponencial

        total += soma

    return total


def integrate_pytorch():

    L = 400*10**(-6)
    Z = np.linspace(0, L, 100)
    qs = np.arange(-75, 75)

    inicio = time.perf_counter()
    
    for q in qs: 
        Aq_pytorch(q, L)

    fim = time.perf_counter()

    total = fim - inicio
    print(f"tempo total integra pytorch: {total} segundos")
    return total


def integrate_scipy():

    L = 400*10**(-6)
    Z = np.linspace(0, L, 100)
    qs = np.arange(-75, 75)

    tempos =[]

    inicio = time.perf_counter()
    
    for q in qs:
        Aq(q, L)
    

    fim = time.perf_counter()

    total = fim - inicio
    print(f"{total}")
    return total

def calcula_tempos():
    count = 0
    tempos = []
    for i in range(0, 50):

        total = integrate_pytorch()
        #total = integrate_scipy()
        tempos.append(total)
    media = (np.sum(tempos))/len(range(0,50))
    print(f"a média de tempo para pytorch foi: {media}, para {len(range(0,50))} execuções")

calcula_tempos()