import os
import time

import numpy as np
import sympy as sym

import scipy as sp 
import scipy.integrate as integrate
import scipy.special as s

import cmath
import math

from mpmath import fac, exp, sin, cos, besselj, legenp, sqrt

import matplotlib.pyplot as plt
plt.rcParams['text.usetex'] = True

import os.path
path = r'D:/Documentos/USP/TCC/Simulações - TCC - BesselBeam/besselBeam/cache/'
import json

from cte import * 
from coefficients import * 
from planeWave import * 


def plot_pw_j1_vs_size():

    x = np.linspace(0.01, 20, 100)

    resultsBlueJ1 = []
    resultsRedJ1 = []
    resultsBlackJ1 = []
    for xi in x:
        #J1(x, M, epislon, mu)
        resultsBlueJ1.append(J1(xi, M, epsilonR2(M, 1), 1))
        resultsRedJ1.append(J1(xi, M2_laop, epsilonR2(M2_laop, 1), 1))
        resultsBlackJ1.append(J1(xi, M3_laop, epsilonR2(M3_laop, 1), 1))

    j1_PW_x_M = {"j1_PW_x_M":[]}
    j1_PW_x_M2_laop = {"j1_PW_x_M2_laop":[]}
    j1_PW_x_M3_laop = {"j1_PW_x_M3_laop":[]}
    
    #gauss bean M = 1.57 - 0.038j
    '''
    for xi in x:

        j1_PW_x_M['j1_PW_x_M'].append({
            "j1": J1(xi, M, epsilonR2(M, 1), 1),
            "x": xi
        })
        json.dump(j1_PW_x_M, open(os.path.abspath(path + 'j1_PW_x_M.json'), 'w'))

        j1_PW_x_M2_laop['j1_PW_x_M2_laop'].append({
            "j1": J1(xi, M2_laop, epsilonR2(M2_laop, 1), 1),
            "x": xi
        })
        json.dump(j1_PW_x_M2_laop, open(os.path.abspath(path + 'j1_PW_x_M2_laop.json'), 'w'))

        j1_PW_x_M3_laop['j1_PW_x_M3_laop'].append({
            "j1": J1(xi, M3_laop, epsilonR2(M3_laop, 1), 1),
            "x": xi
        })
        json.dump(j1_PW_x_M3_laop, open(os.path.abspath(path + 'j1_PW_x_M3_laop.json'), 'w'))
        
    '''

    #plt.figure(figsize=[7,5])
    plt.plot(x, resultsBlueJ1, 'b--', label= "OP ($g_n$ = 1)")
    plt.plot(x, resultsRedJ1, 'r-.', label= "OP")
    plt.plot(x, resultsBlackJ1, 'g', label= "OP")
    plt.xlabel('Size Parameter x')
    plt.ylabel('Asymmetry Factor J1')
    plt.grid()
    plt.legend(loc='best')
    plt.text(12.5, .03, 'M = 1.57 - 0.038j')
    plt.text(5, -.1, 'M = 1.57 - 0.19j')
    plt.text(8, -.33, 'M = 1.57 - 0.95j')
    plt.savefig(r'D:/Documentos/USP/TCC/Simulações - TCC - BesselBeam/besselBeam/graphs/PW_j1_vs_size.png', format='png', dpi=500)
    plt.show()

plot_pw_j1_vs_size()