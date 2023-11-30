import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from mpmath import exp
from mpl_toolkits.mplot3d import Axes3D
import scipy.special as s
plt.rcParams['text.usetex'] = True

lamb =1064*10**(-9)
#lamb = 9.58*10**(-6)
k = (2 * np.pi)/lamb
L = 400*10**(-6)
Q = 0.8*k
M = (1.57 - 0.038j)
mi = 1
eta_r = 1/M
c = 3.0 * (10**8)
nano  = 10**(-9)
micro = 10**(-6)
mili  = 10**(-3)
cm = 10**(-2)
deltaRho0 = 60 * micro
l = 1064 * nano
c = 3.0 * (10**8)
k0 = 2*math.pi / l
t = 1

def k_rho(k, alpha):
    return (k*np.sin(alpha))

def k_z(k, alpha):
    return (k*np.cos(alpha))

spot5 = ((405/100)/(k * np.sin(math.radians(5))))
spot10 = ((405/100)/(k * np.sin(math.radians(10))))
spot15 = ((405/100)/(k * np.sin(math.radians(15))))

def psiBessel(rho, alpha, z_0):
    i = 1j
    cosAxicon = np.cos(alpha)
    bessel = (s.jv(0, k*np.sin(alpha)*rho))
    expo = (np.exp(-1j*k*cosAxicon*z_0))

    resul = bessel

    return resul

'''
# def psiGauss(rho, z, t, a, k0, c):
qntPoints = 200
z = np.linspace(0, 1.25, 200)
results = np.zeros((qntPoints,qntPoints))

rho = np.linspace(100, 100, 200)
phi = np.linspace(0, 2*np.pi, 200)
R, P = np.meshgrid(rho, phi)

for rho_i in range(qntPoints):
    for zi in range(qntPoints):
        results[rho_i][zi]= (psiBessel(rho[rho_i], math.radians(5), z[zi]))
        #resultsZ09.append((psiBessel(i, math.radians(10), zi)))
        #resultsZ18.append((psiBessel(i, math.radians(15), zi)))


# Express the mesh in the cartesian system.
X, Y = R*np.cos(P), R*np.sin(P)

'''
#pToPrint = np.linspace(-0.4, 0.4, 200)
#PToPrint, ZToPrint = np.meshgrid(pToPrint, z, indexing='ij')
'''

fig = plt.figure(figsize=(12, 7))
fig3d = fig.add_subplot(projection='3d')

fig3d.plot_surface(X, Y, results)
#plt.plot(p, resultsZ0, 'b', label=r'$\alpha = 5^{\circ}$')
#plt.plot(p, resultsZ09, 'g--', label=r'$\alpha = 10^{\circ}$')
#plt.plot(p, resultsZ18, 'r.', label=r'$\alpha = 15^{\circ}$')

plt.title(r'Feixe de Bessel para $\alpha$ arbitrários')
#plt.xlabel(r'$\rho$ (mm)')
plt.ylabel(r'$z$ (mm)')
plt.grid()
plt.legend()
plt.show()

#plt.savefig(r'D:/Documentos/USP/TCC/Simulações - TCC - BesselBeam/besselBeam/graphs/psi_fg_p.png', format='png', dpi=1000)
'''

def plot_bessel_func():

    X = np.linspace(0, 25, 200)
    bessel0 = []
    bessel1 = []
    bessel2 = []
    bessel3 = []
    bessel10 = []
    bessel50 = []

    for xi in X:
        bessel0.append(s.j0(xi))
        bessel1.append(s.j1(xi))
        bessel2.append(s.jn(5, xi))
        #bessel3.append(s.jn(3, xi))
        bessel10.append(s.jn(10, xi))
        bessel50.append(s.jn(15, xi))

    
    plt.figure(figsize=(12, 5))
    plt.plot(X, bessel0, 'y', label=r'$\mathcal{J}_{0}$')
    #plt.plot(X, bessel2, 'b', label=r'$\mathcal{J}_{2}$')
    plt.plot(X, bessel2, 'm', label=r'$\mathcal{J}_{5}$')
    plt.plot(X, bessel10, 'c', label=r'$\mathcal{J}_{10}$')
    plt.plot(X, bessel50 , 'k-.', label=r'$\mathcal{J}_{15}$')
    #plt.title(r'Função de Bessel para $n$ inteiros')
#plt.xlabel(r'$\rho$ (mm)')
    plt.ylabel(r'$\mathcal{J}(x)$')
    plt.grid()
    plt.legend()
    plt.show()


plot_bessel_func()