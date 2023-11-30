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

from legendre import *
from specials import * 
from besselBeam.cte import * 
from Aq import * 
from beamShape import *
from coefficients import * 
from planeWave import * 


'''
print('Versão Serial')
start_time = time.time()
A_mn( 2, M, 0.1, 0, pi/20, 0.1*(((405/100)/(k * math.sin(pi/20)))), 0, 0, 0, 0)
end_time = time.time()
print('Tempo de execução (versão serial): {} segundos'.format(end_time - start_time))
'''
ceillingX = lambda x : math.ceil(x + (4.05 * (x ** (1/3))) + 2)
print(ceillingX(3))

