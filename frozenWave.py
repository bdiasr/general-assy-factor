import os
import time

import numpy as np

import scipy as sp 
import scipy.integrate as integrate
import scipy.special as s

import cmath
import math


import torch
from torchquad import Simpson, set_up_backend
import torch.nn as nn


import matplotlib.pyplot as plt
plt.rcParams['text.usetex'] = True
import matplotlib as mpl

from cte import *
from coefficients import * 

from joblib import Parallel, delayed

import warnings
warnings.filterwarnings(action='ignore', category=UserWarning)



device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
set_up_backend("torch", data_type="float64")
num_gpus = torch.cuda.device_count()
torch.set_printoptions(precision=20)
dimension = 1
simp = Simpson()
integration_domain = [[0, L]] * dimension
Qs = np.arange(-75, 76)
points = 130


def F_z1(z, L):
    z_tensor = torch.as_tensor(z)
    L_tensor = torch.as_tensor(L)
    
    condition = (z_tensor >= 3*L_tensor/8) & (z_tensor <= 5*L_tensor/8)

    value = torch.where(condition, 
                        torch.tensor(1),
                        torch.tensor(0))
    return(value)


def F_z2(z, L):

    z_tensor = torch.as_tensor(z)
    L_tensor = torch.as_tensor(L)

    condition = ((z_tensor >= 2*L_tensor/8) & ( z_tensor <= 6*L_tensor/8))
    
    value = torch.where(condition,
                        torch.exp(((10**(4))*(z_tensor - 6*L_tensor/8)/2)),
                        torch.tensor(0))
    
    return(value)


def F_z3(z, L):

    z_tensor = torch.as_tensor(z)
    L_tensor = torch.as_tensor(L)

    condition = (z_tensor >= L_tensor / 12) & (z_tensor <= 11 * L_tensor / 12)
    expoente = -5 * ((z_tensor - (0.5 * L_tensor)) ** 2 * 1 / L_tensor ** 2)

    value = torch.where(condition,
                        torch.exp(expoente) * torch.cos(6 * torch.pi * z_tensor / L_tensor),
                        torch.tensor(0))
    return value

def F_z4(z, L):
    z_tensor = torch.as_tensor(z)
    L_tensor = torch.as_tensor(L)

    condition = ((z_tensor >= 1/4*L_tensor) & ( z_tensor <= 3/4*L_tensor))
    
    value = torch.where(condition,
                        torch.tensor((torch.special.bessel_j0(1.6*10**(-6) * z_tensor))),
                        torch.tensor(0))
    
    return(value)
    


def complex_integrate(L, q, a, b):

    integration_domain = [[a, b]]
    
    integrand = lambda z: (F_z2(z, L) * 
                           torch.exp(torch.as_tensor(1j*z*(2*torch.pi*q)/L)))
    
    real_integrand = lambda x: integrand(x).real
    complex_integrand = lambda x: integrand(x).imag
    
    real_integral = simp.integrate((real_integrand), dim=1, N=301, integration_domain=integration_domain).item()
    complex_integral = simp.integrate((complex_integrand), dim=1, N=301, integration_domain=integration_domain).item()
    
    if complex_integral == 0.0:
        integral_val = real_integral
    else:
        integral_val = (real_integral + 1j*complex_integral)
        
    return integral_val


def Aq(q, L):

    frac = 1/L
    integ = complex_integrate(L, q, 0, L)

    return (frac*integ)

def Psi_args(rho, z, q):

    Q = 0.8*k
    L = 400*10**(-6)
    soma = []

    k_zq = Q + 2*np.pi*q/L
    k_pq = np.sqrt(k**2 - k_zq**2)
    a = Aq(q, L)
    j0 = (s.j0((k_pq *rho)))
    exponencial = np.exp((-1j * k_zq * z))
    soma = a * j0 * exponencial
    
    #total += soma
    return soma

def Psi_Parallel(rho, z):

    soma = []
    total = 0
    qs = np.arange(-75, 76)
    soma.append(Parallel(n_jobs=1)(delayed(Psi_args)(rho, z, q) for q in qs))
    total = np.sum(soma)
    return total

def Psi(rho, z):
    Q = 0.8*k
    L = 400*10**(-6)
    soma = []
    total = 0
    qs = np.arange(-75, 76)
    
    for q in qs:

        k_zq = Q + 2*np.pi*q/L
        k_pq = np.sqrt(k**2 - k_zq**2)
        a = Aq(q, L)
        j0 = (s.j0((k_pq *rho)))
        exponencial = np.exp((-1j * k_zq * z))
        soma = a * j0 * exponencial
      
        total += soma

    return total


def pi_n(m, n, x):
    
    fator = cmath.sqrt(1 - x**2)
    frac =  1/fator
    aux = s.lpmn(m, n, x)[0]
    pi_val = aux[m][n]*frac
    
    return pi_val

def tau_n(m, n, x):

    fator = cmath.sqrt(1 - x**2)
    aux = s.lpmn(m, n, x)[1]
    tau_val = -(fator)*aux[m][n]
    return tau_val

def gn_FrozenWave(n, N, k, L, Q, z0):

    soma = []
    total = 0

    primeiroTermoDem = n*(n + 1)
    qs = np.arange(N, (-N + 1))

    for q in qs:

        k_zq = Q + ((2*np.pi*q)/L)
        k_termo = k_zq/k
        primeiroTermoSum = Aq(q, L)/(1 + k_termo)
        primeiroTermoMul = (pi_n(1, n, k_termo)) + (tau_n(1, n, k_termo))
        exponencial = np.exp(1j * k_zq * z0)

        soma = primeiroTermoSum * primeiroTermoMul * exponencial
        total += soma
  
    gn = (-2/primeiroTermoDem)*(total)

    return gn

n_max = (lambda x: math.ceil(x + 4.05*x**(1/3)+2))


def j_FW_args(epslon2, m, x, z0, mi, k, i):

    N = -75
    L = 400*10**(-6)
    Q = 0.8*k

    soma=[]

    cn = cs_n(m, mi, i, x)
    conjCn = np.conjugate(cn)
    cn1 = cs_n(m, mi, i+1, x)
    conjCn1 = np.conjugate(cn1)
    gnVar = gn_FrozenWave(i, N, k, L, Q, z0)
    conjGn = np.conjugate(gnVar)
    gn1Var = gn_FrozenWave(i+1, N, k, L, Q, z0)
    conjGn1 = np.conjugate(gn1Var)
    
    rn = r_n(m, i, x)
    rn1 = r_n(m, i+1, x)
    
    dn = ds_n(m, mi, i, x)
    conjDn = np.conjugate(dn)
    dn1 = ds_n(m, mi, i+1, x)
    sn = S_n(m, i, x)
    b = (i*(i+2)/m)
    c1 = (gn1Var)*(conjGn)*(cn1)*(conjCn)*(rn1)
    c2 = (((abs(1/m)**2))*(gn1Var)*(conjGn)*(dn1)*(conjDn)*(rn))
    
    d = (i*(i+2)/(i+1))
    d1 = (gnVar)*(conjGn1)*(cn)*(conjCn1)
    d2 = ((abs(1/m)**2)*(gn1Var)*(conjGn)*(dn1)*(conjDn))
    f = ((2*i+1)/(i*(i+1))) * (abs(gnVar)**2)*(cn*conjDn)

    soma = b*(c1+c2) - (((d*(d1 + d2))+ f*(1/m))*S_n(m, i, x))

    return soma

def j_FW_Parallel(epslon2, m, x, z0, mi, k):

    n_maximo = n_max(x)
    a = (3*epslon2)/((abs(m)**2)*(x**3))
    soma = []
    total = 0

    soma.append(Parallel(n_jobs=2)(delayed(j_FW_args)(epslon2, m, x, z0, mi, k, i) for i in range(1, n_maximo + 1)))
    total = np.sum(soma)
    j = a * total.imag
    #print("j:", j, "z:", z0)
    return j

def j_any(epslon2, m, x, z0, mi, k):
    
    n_maximo = n_max(x)
    
    a = (3*epslon2)/((abs(m)**2)*(x**3))
    soma = []
    total = 0
    N = -75
    L = 400*10**(-6)
    Q = 0.8*k


    for i in range(1, (n_maximo + 1)):
        
        cn = cs_n(m, mi, i, x)
        conjCn = np.conjugate(cn)
        
        cn1 = cs_n(m, mi, i+1, x)
        conjCn1 = np.conjugate(cn1)
        
        gnVar = gn_FrozenWave(i, N, k, L, Q, z0)
        conjGn = np.conjugate(gnVar)
        
        gn1Var = gn_FrozenWave(i+1, N, k, L, Q, z0)
        conjGn1 = np.conjugate(gn1Var)
        
        rn = r_n(m, i, x)
        rn1 = r_n(m, i+1, x)
        
        dn = ds_n(m, mi, i, x)
        conjDn = np.conjugate(dn)
        
        dn1 = ds_n(m, mi, i+1, x)
        sn = S_n(m, i, x)
        
        b = (i*(i+2)/m)
        c1 = (gn1Var)*(conjGn)*(cn1)*(conjCn)*(rn1)
        
        c2 = (((abs(1/m)**2))*(gn1Var)*(conjGn)*(dn1)*(conjDn)*(rn))
        
        d = (i*(i+2)/(i+1))
        d1 = (gnVar)*(conjGn1)*(cn)*(conjCn1)
        d2 = ((abs(1/m)**2)*(gn1Var)*(conjGn)*(dn1)*(conjDn))
       
        f = ((2*i+1)/(i*(i+1))) * (abs(gnVar)**2)*(cn*conjDn)

        soma = b*(c1+c2) - (((d*(d1 + d2))+ f*(1/m))*S_n(m, i, x))
        
        total += soma
    j = a * total.imag
        
    return j


'''
# usar como referencia 
epslon1 = (M**2).real
epslon2 = -(M**2).imag

epslon22 = -(M2**2).imag
epslon23 = -(M3**2).imag
'''


def plot_Psi():

    """ Figure DPI """
    TEST_DPI = 100
    FINISH_DPI = 200

    mpl.rcParams["lines.linewidth"] = 2
    mpl.rcParams["text.usetex"] = True
    mpl.rcParams["font.family"] = "serif"
    mpl.rcParams["axes.titlesize"] = 32
    mpl.rcParams["axes.labelsize"] = 24
    mpl.rcParams["axes.titlepad"] = 24
    mpl.rcParams["xtick.labelsize"] = "xx-large"
    mpl.rcParams["ytick.labelsize"] = "xx-large"
    mpl.rcParams["legend.fontsize"] = 24
    mpl.rcParams['figure.figsize'] = [12, 8]
    mpl.rcParams['figure.dpi'] = TEST_DPI

    Psi_values = []
    L = 400*10**(-6)

    Z = np.linspace(0, L, 100)

    print('Plot de Psi com Pytorch - sequencial')
    start_time = time.time()

    for z in Z: 
        Psi_values.append(abs(Psi(0, z))**2)

    end_time = time.time()
    print('Tempo de execução (versão pytorch - sequencial): {} segundos'.format(end_time - start_time))

    plt.plot(Z, Psi_values, 'b')
    plt.title(r'$|\Psi (0, z)|^2$')
    plt.xlabel(r'$z (\mu m)$')
    plt.grid()
    plt.show()


def plot_Psi_Parallel():

    """ Figure DPI """
    TEST_DPI = 100
    FINISH_DPI = 200

    mpl.rcParams["lines.linewidth"] = 2
    mpl.rcParams["text.usetex"] = True
    mpl.rcParams["font.family"] = "serif"
    mpl.rcParams["axes.titlesize"] = 32
    mpl.rcParams["axes.labelsize"] = 24
    mpl.rcParams["axes.titlepad"] = 24
    mpl.rcParams["xtick.labelsize"] = "xx-large"
    mpl.rcParams["ytick.labelsize"] = "xx-large"
    mpl.rcParams["legend.fontsize"] = 24
    mpl.rcParams['figure.figsize'] = [12, 8]
    mpl.rcParams['figure.dpi'] = TEST_DPI

    Psi_values = []
    L = 400*10**(-6)

    Z = np.linspace(0, L, 100)

    print('Plot de Psi com Pytorch - paralela')
    start_time = time.time()

    for z in Z: 
        Psi_values.append(abs(Psi_Parallel(0, z))**2)

    end_time = time.time()
    print('Tempo de execução (versão pytorch - paralela): {} segundos'.format(end_time - start_time))

    plt.plot(Z, Psi_values, 'b')
    plt.title(r'$|\Psi (0, z)|^2$')
    plt.xlabel(r'$z (\mu m)$')
    plt.grid()
    plt.show()

def plot_Fz():
    
    """ Figure DPI """
    TEST_DPI = 100
    mpl.rcParams["lines.linewidth"] = 2
    mpl.rcParams["text.usetex"] = True
    mpl.rcParams["font.family"] = "serif"
    mpl.rcParams["axes.titlesize"] = 32
    mpl.rcParams["axes.labelsize"] = 24
    mpl.rcParams["axes.titlepad"] = 24
    mpl.rcParams["xtick.labelsize"] = "xx-large"
    mpl.rcParams["ytick.labelsize"] = "xx-large"
    mpl.rcParams["legend.fontsize"] = 24
    mpl.rcParams['figure.figsize'] = [12, 8]
    mpl.rcParams['figure.dpi'] = TEST_DPI # Use 200 for fine images - note this is slow
    #mpl.rcParams["axes.grid"] = True

    
    F_values = []
    L = (400*10**(-6))

    Z = torch.linspace(0, L, 100)

    print('Plot de F_z com Pytorch')
    start_time = time.time()

    for z in Z: 
        F_values.append(abs(F_z1(z, L))**2)

    end_time = time.time()
    print('Tempo de execução (versão pytorch): {} segundos'.format(end_time - start_time))

    plt.plot(Z.cpu(), torch.as_tensor(F_values).cpu(), 'r')
    plt.xlabel(r'$z$')
    plt.title(r'$F(z)$')
    plt.grid()
    plt.show()


def plot_Jz_vs_z():
    
    zval = np.linspace(0, 400*10**(-6), 15)
    z_axis = np.linspace(0, 400, 15)

    cnt = 0
    jn1_gauss_blue = []
    jn2_gauss_red = []
    jn3_gauss_green = []
   
    start_time = time.time()
    #M = 1.57 - 0.038j
    for i in zval:
        print("i:", cnt)
        jn1_gauss_blue.append(j_any(epsilonR2(M, mi), M, 0.1, i, mi, k)*250)
        #jn2_gauss_red.append(j_any(epsilonR2(M, mi), M, 3, i, mi, k))
        #jn3_gauss_green.append(j_any(epsilonR2(M, mi), M, 8, i, mi, k))
        cnt = cnt + 1
    end_time = time.time()
    print('Plot de J1 vs z_0 - Tempo de execução (versão sequencial): {} segundos'.format(end_time - start_time))

    #plt.figure(figsize=[10,8])
    plt.plot(z_axis, jn1_gauss_blue,"b", label= "x = 0.1" )
    #plt.plot(z_axis, jn2_gauss_red,"r", label= "x = 3" )    
    #plt.plot(z_axis, jn3_gauss_green,"g", label= "x = 8" )

    plt.xlabel(r'$z_0$ $\mu m$')
    plt.ylabel(r'Assymetry Factor $J_1$')

    plt.legend()
    plt.grid()
    plt.show()

def plot_J1_FW_z_Parallel():
    
    zval = np.linspace(0, 400*10**(-6), 15)
    z_axis = np.linspace(0, 400, 15)

    cnt = 0
    jn1_gauss_blue = []
    jn2_gauss_red = []
    jn3_gauss_green = []

    start_time = time.time()
   
    #M = 1.57 - 0.038j
    for i in zval:
        print("i:", cnt)
        jn1_gauss_blue.append(j_FW_Parallel(epsilonR2(M, mi), M, 0.1, i, mi, k)*250)
        #jn2_gauss_red.append(j_FW_Parallel(epsilonR2(M, mi), M, 3, i, mi, k))
        #jn3_gauss_green.append(j_FW_Parallel(epsilonR2(M, mi), M, 8, i, mi, k))
        cnt = cnt + 1
    end_time = time.time()
    print('Plot de J1 vs z_0 - Tempo de execução (versão paralela): {} segundos'.format(end_time - start_time))

    #plt.figure(figsize=[10,8])
    plt.plot(z_axis, jn1_gauss_blue,"b", label= "x = 0.1" )
    #plt.plot(z_axis, jn2_gauss_red,"r", label= "x = 3" )    
    #plt.plot(z_axis, jn3_gauss_green,"g", label= "x = 8" )

    plt.xlabel(r'$z_0$ $\mu m$')
    plt.ylabel(r'Assymetry Factor $J_1$')

    plt.legend()
    plt.grid()
    plt.show()

'''
plot_Fz()
plot_Psi()
'''
#plot_Psi()
#plot_Psi_Parallel()
#plot_J1_FW_z_Parallel()
#plot_Jz_vs_z()
