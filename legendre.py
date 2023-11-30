import numpy as np
import os
import mpmath
from mpmath import fac, exp, sin, cos, besselj, legenp, sqrt
import time


import matplotlib.pyplot as plt
plt.rcParams['text.usetex'] = True


from mpmath import mp, fp

# Definindo a precisão da aritmética de ponto flutuante
fp.prec = 30 # Número de bits significativos em ponto flutuante
# Definindo a precisão decimal
mp.dps = 80 # Número de dígitos decimais

def pi_mn(m, n, x):
    
    cosAlpha = cos(x)
    sinAlpha = sin(x)
    
    num = legenp(n, m, cosAlpha)
    if(sinAlpha == 0):
        return 0
    else:
        return num/sinAlpha

def tau_mn(m, n, x):
    
    #x = cos de alpha 

    pmn = legenp(n, m, x)
    pmn1 =  legenp(n+1, m, x)

    fristTerm = -(n + 1)*x*pmn
    secondTerm = (n - m + 1)*pmn1

    num = fristTerm + secondTerm
    
    if x == 1:
        den = sqrt(1 - 0.999999999)
    else:
        den = sqrt(1 - (x**2))
    
    return num/den


'''
print('Versão Mpmath')
start_time = time.time()
print(pi_mn(1, 1, pi/360))
print(tau_mn(1, 1, cos(pi/360)))
end_time = time.time()
print('Tempo de execução (Mpmath): {:.6f} segundos'.format(end_time - start_time))
'''

def fig_tau():

    n1 = []
    n2 = []
    n3 = []

    X = np.linspace(-0.99, 0.99, 200)

    for x in X:
        n1.append(tau_mn(1, 1, cos(x)))

    for x in X:
        n2.append(tau_mn(1, 2, cos(x)))
  
    for x in X:
        n3.append(tau_mn(1, 3, cos(x)))

    plt.figure(figsize=[7,5])
    plt.plot(X, n1, 'b', label=r'$n_1$')
    plt.plot(X, n2, 'r', label=r'$n_3$')
    plt.plot(X, n3, 'g', label=r'$n_3$')

    plt.ylim(-6, 6)
    plt.xlabel(r'$\cos{\alpha}$')
    #plt.ylabel('')
    plt.title(r'$\tau^{m}_{n}(\cos{\alpha})$')
    plt.legend()
    plt.grid()
    plt.show()


def fig_pi():

    n1 = []
    n2 = []
    n3 = []

    X = np.linspace(-1, 1, 200)

    for x in X:
        n1.append(pi_mn(1, 1, x))

    for x in X:
        n2.append(pi_mn(1, 2, x))
  
    for x in X:
        n3.append(pi_mn(1, 3, x))

    plt.plot(X, n1, 'b', label ='n1')
    plt.plot(X, n2, 'r', label='n2')
    plt.plot(X, n3, 'g', label='n3')

    #plt.ylim(-6, 2)
    plt.xlabel('x')
    #plt.ylabel('')
    plt.legend()
    plt.grid()
    plt.show()

#fig_tau()
#fig_pi()