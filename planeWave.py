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
from cte import * 
from Aq import * 
from beamShape import * 
from coefficients import * 

ceillingX = lambda x : math.ceil(x + (4.05 * (x ** (1/3))) + 2)

def J1(x, M, epislon, mu):

    firstTerm = (3 * epislon) / ((abs(M) ** 2) * (x ** 3))
    summationResult = 0
    nMax = ceillingX(x)
    etaR = 1 / M
    mod2EtaR = abs(etaR) ** 2
    
    for n in range(1, (nMax+1)):
        cn = cs_n(M, mi, n, x) 
        conjCn = np.conj(cn)
        cn1 = cs_n(M, mi, n + 1, x) 
        conjCn1 = np.conj(cn1)

        rn = r_n(M, n, x)
        rn1 = r_n(M, n + 1, x)

        dn = ds_n(M, mi, n, x)
        conjDn = np.conj(dn)
        dn1 =  ds_n(M, mi, n + 1, x)
        conjDn1 = np.conj(dn1)
        sn = S_n(M, n, x)
        
        fristTermFor = ((n*(n+2)) / M) * ((cn1*conjCn*rn1) + (mod2EtaR*dn1*conjDn*rn))
        secondTermFor = ((n*(n+2)) / (n+1)) * ((cn*conjCn1) + (mod2EtaR*dn1*conjDn))
        thirdTermFor =  np.conj(etaR) * (((2*n)+1) / (n*(n+1))) * (cn*conjDn)
        
        result = fristTermFor - ((secondTermFor + thirdTermFor) * sn)
        
        summationResult += result
        
    return firstTerm * summationResult.imag