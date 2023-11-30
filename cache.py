import os
import time
import json 

import numpy as np



import matplotlib.pyplot as plt
plt.rcParams['text.usetex'] = True
import matplotlib as mpl
mpl.rcParams['text.usetex'] = True
mpl.rcParams['text.latex.preamble'] = r'\usepackage{{amsmath}}'

from cte import *
from coefficients import * 
from frozenWave import * 
from asymetryVector import *
from planeWave import * 

import os.path
path = r'D:/Documentos/USP/TCC/Simulações - TCC - BesselBeam/besselBeam/cache/'


def cache_FW1_Psi():
    
    zval = np.linspace(0, 400*10**(-6), 200)

    cnt = 0
    FW1_Psi = {"FW1_Psi":[]}
    #adicionar verificação se o arquivo ja existe 

    #M = 1.57 - 0.038j
    for zi in zval:
        print("i:", cnt)
        FW1_Psi['FW1_Psi'].append({
            "Psi": abs(Psi(0, zi))**2,
            "z": zi
        })
        json.dump(FW1_Psi, open(os.path.abspath(path + 'FW1_Psi.json'), 'w'))
        cnt = cnt + 1

def cache_FW2_Psi():
    
    zval = np.linspace(0, 400*10**(-6), 200)

    cnt = 0
    FW2_Psi = {"FW2_Psi":[]}
    #adicionar verificação se o arquivo ja existe 

    #M = 1.57 - 0.038j
    for zi in zval:
        print("i:", cnt)
        FW2_Psi['FW2_Psi'].append({
            "Psi": abs(Psi(0, zi))**2,
            "z": zi
        })
        json.dump(FW2_Psi, open(os.path.abspath(path + 'FW2_Psi.json'), 'w'))
        cnt = cnt + 1


def cache_FW3_Psi():
    
    zval = np.linspace(0, 400*10**(-6), 200)

    cnt = 0
    FW3_Psi = {"FW3_Psi":[]}
    #adicionar verificação se o arquivo ja existe 

    #M = 1.57 - 0.038j
    for zi in zval:
        print("i:", cnt)
        FW3_Psi['FW3_Psi'].append({
            "Psi": abs(Psi(0, zi))**2,
            "z": zi
        })
        json.dump(FW3_Psi, open(os.path.abspath(path + 'FW3_Psi.json'), 'w'))
        cnt = cnt + 1

def cache_FW4_Psi():
    
    zval = np.linspace(0, 400*10**(-6), 200)

    cnt = 0
    FW3_Psi = {"FW4_Psi":[]}
    #adicionar verificação se o arquivo ja existe 

    #M = 1.57 - 0.038j
    for zi in zval:
        print("i:", cnt)
        FW3_Psi['FW4_Psi'].append({
            "Psi": abs(Psi(0, zi))**2,
            "z": zi
        })
        json.dump(FW3_Psi, open(os.path.abspath(path + 'FW4_Psi.json'), 'w'))
        cnt = cnt + 1
 
def cache_FW1_J1_z():
    
    zval = np.linspace(0, 400*10**(-6), 200)

    cnt = 0
    j1_FW_01 = {}
    j1_FW_3 = {}
    j1_FW_8 = {}

    j1_FW_01 = {"j1_FW_01":[]}
    j1_FW_3 = {"j1_FW_3":[]}
    j1_FW_8 = {"j1_FW_8":[]}

    #adicionar verificação se o arquivo ja existe 

    #M = 1.57 - 0.038j
    for zi in zval:
        print("i:", cnt)

        j1_FW_01['j1_FW_01'].append({
            "j1": j_any(epsilonR2(M, mi), M, 0.1, zi, mi, k),
            "z": zi
        })
        json.dump(j1_FW_01, open(os.path.abspath(path + 'FW_J1_01_z.json'), 'w'))

        j1_FW_3['j1_FW_3'].append({
            "j1": j_any(epsilonR2(M, mi), M, 3, zi, mi, k),
            "z": zi
        })
        json.dump(j1_FW_3, open(os.path.abspath(path + 'FW_J1_3_z.json'), 'w'))

        j1_FW_8['j1_FW_8'].append({
            "j1": j_any(epsilonR2(M, mi), M, 8, zi, mi, k),
            "z": zi
        })
        json.dump(j1_FW_8, open(os.path.abspath(path + 'FW_J1_8_z.json'), 'w'))

        cnt = cnt + 1

def cache_FW2_J1_z():
    
    zval = np.linspace(0, 400*10**(-6), 200)
    cnt = 0

    j1_FW2_01 = {"j1_FW2_01":[]}
    j1_FW2_3 = {"j1_FW2_3":[]}
    j1_FW2_8 = {"j1_FW2_8":[]}

    #adicionar verificação se o arquivo ja existe 

    #M = 1.57 - 0.038j
    for zi in zval:
        print("i:", cnt)

        j1_FW2_01['j1_FW2_01'].append({
            "j1": j_any(epsilonR2(M, mi), M, 0.1, zi, mi, k),
            "z": zi
        })
        json.dump(j1_FW2_01, open(os.path.abspath(path + 'FW2_J1_01_z.json'), 'w'))

        j1_FW2_3['j1_FW2_3'].append({
            "j1": j_any(epsilonR2(M, mi), M, 3, zi, mi, k),
            "z": zi
        })
        json.dump(j1_FW2_3, open(os.path.abspath(path + 'FW2_J1_3_z.json'), 'w'))

        j1_FW2_8['j1_FW2_8'].append({
            "j1": j_any(epsilonR2(M, mi), M, 8, zi, mi, k),
            "z": zi
        })
        json.dump(j1_FW2_8, open(os.path.abspath(path + 'FW2_J1_8_z.json'), 'w'))

        cnt = cnt + 1

def cache_FW3_J1_z():
    
    zval = np.linspace(0, 400*10**(-6), 200)
    cnt = 0

    j1_FW3_01 = {"j1_FW3_01":[]}
    j1_FW3_3 = {"j1_FW3_3":[]}
    j1_FW3_8 = {"j1_FW3_8":[]}

    #adicionar verificação se o arquivo ja existe 

    #M = 1.57 - 0.038j
    for zi in zval:
        print("i:", cnt)

        j1_FW3_01['j1_FW3_01'].append({
            "j1": j_any(epsilonR2(M, mi), M, 0.1, zi, mi, k),
            "z": zi
        })
        json.dump(j1_FW3_01, open(os.path.abspath(path + 'FW3_J1_01_z.json'), 'w'))

        j1_FW3_3['j1_FW3_3'].append({
            "j1": j_any(epsilonR2(M, mi), M, 3, zi, mi, k),
            "z": zi
        })
        json.dump(j1_FW3_3, open(os.path.abspath(path + 'FW3_J1_3_z.json'), 'w'))

        j1_FW3_8['j1_FW3_8'].append({
            "j1": j_any(epsilonR2(M, mi), M, 8, zi, mi, k),
            "z": zi
        })
        json.dump(j1_FW3_8, open(os.path.abspath(path + 'FW3_J1_8_z.json'), 'w'))

        cnt = cnt + 1

def cache_FW4_J1_z():
    
    zval = np.linspace(0, 400*10**(-6), 200)
    cnt = 0

    j1_FW3_01 = {"j1_FW4_01":[]}
    j1_FW3_3 = {"j1_FW4_3":[]}
    j1_FW3_8 = {"j1_FW4_8":[]}

    #adicionar verificação se o arquivo ja existe 

    #M = 1.57 - 0.038j
    for zi in zval:
        print("i:", cnt)

        j1_FW3_01['j1_FW4_01'].append({
            "j1": j_any(epsilonR2(M, mi), M, 0.1, zi, mi, k),
            "z": zi
        })
        json.dump(j1_FW3_01, open(os.path.abspath(path + 'FW4_J1_01_z.json'), 'w'))

        j1_FW3_3['j1_FW4_3'].append({
            "j1": j_any(epsilonR2(M, mi), M, 3, zi, mi, k),
            "z": zi
        })
        json.dump(j1_FW3_3, open(os.path.abspath(path + 'FW4_J1_3_z.json'), 'w'))

        j1_FW3_8['j1_FW4_8'].append({
            "j1": j_any(epsilonR2(M, mi), M, 8, zi, mi, k),
            "z": zi
        })
        json.dump(j1_FW3_8, open(os.path.abspath(path + 'FW4_J1_8_z.json'), 'w'))

        cnt = cnt + 1
    
def cache_FW1_J1_x():

    xval = np.linspace(0.1, 20, 100)
    cnt = 0

    j1_FW1_x_L2 = {"j1_FW1_x_L2":[]}
    j1_FW1_x_L4 = {"j1_FW1_x_L4":[]}
    
    #gauss bean M = 1.57 - 0.038j
    for xi in xval:
        print("i:", cnt)

        j1_FW1_x_L2['j1_FW1_x_L2'].append({
            "j1": j_any(epsilonR2(M, mi), M, xi, L/2, mi, k),
            "x": xi
        })
        json.dump(j1_FW1_x_L2, open(os.path.abspath(path + 'FW1_j1_x_L2.json'), 'w'))
        
        j1_FW1_x_L4['j1_FW1_x_L4'].append({
            "j1": j_any(epsilonR2(M, mi), M, xi, L/4, mi, k),
            "x": xi
        })
        json.dump(j1_FW1_x_L4, open(os.path.abspath(path + 'FW1_j1_x_L4.json'), 'w'))

        cnt = cnt + 1

def cache_FW3_J1_x():

    xval = np.linspace(0.1, 20, 100)
    cnt = 0

    j1_FW3_x_L2 = {"j1_FW3_x_L2":[]}
    j1_FW3_x_L4 = {"j1_FW3_x_L4":[]}
    
    #gauss bean M = 1.57 - 0.038j
    for xi in xval:
        print("i:", cnt)

        j1_FW3_x_L2['j1_FW3_x_L2'].append({
            "j1": j_any(epsilonR2(M, mi), M, xi, L/2, mi, k),
            "x": xi
        })
        json.dump(j1_FW3_x_L2, open(os.path.abspath(path + 'FW3_j1_x_L2.json'), 'w'))
        
        j1_FW3_x_L4['j1_FW3_x_L4'].append({
            "j1": j_any(epsilonR2(M, mi), M, xi, L/4, mi, k),
            "x": xi
        })
        json.dump(j1_FW3_x_L4, open(os.path.abspath(path + 'FW3_j1_x_L4.json'), 'w'))

        cnt = cnt + 1

def cache_FW2_J1_x():

    xval = np.linspace(0.1, 20, 100)
    cnt = 0

    j1_FW2_x_L2 = {"j1_FW2_x_L2":[]}
    j1_FW2_x_L4 = {"j1_FW2_x_L4":[]}
    
    #gauss bean M = 1.57 - 0.038j
    for xi in xval:
        print("i:", cnt)

        j1_FW2_x_L2['j1_FW2_x_L2'].append({
            "j1": j_any(epsilonR2(M, mi), M, xi, L/2, mi, k),
            "x": xi
        })
        json.dump(j1_FW2_x_L2, open(os.path.abspath(path + 'j1_FW2_x_L2.json'), 'w'))
        
        j1_FW2_x_L4['j1_FW2_x_L4'].append({
            "j1": j_any(epsilonR2(M, mi), M, xi, L/4, mi, k),
            "x": xi
        })
        json.dump(j1_FW2_x_L4, open(os.path.abspath(path + 'j1_FW2_x_L4.json'), 'w'))

        cnt = cnt + 1

def cache_BB_J1_x_alphas():

    xval = np.linspace(0.1, 20, 100)
    cnt = 0

    '''
    configuracoes: onda plana, alpha = 1, alpha = 30 e alpha = 60
    '''

    j1_PW_x = {"M": '1.57 - 0.038j', "j1_PW_x":[]}
    j1_BB_x_alpha1 = {"M": "1.57 - 0.038j", "j1_BB_x_alpha1":[]}
    j1_BB_x_alpha30 = {"M": "1.57 - 0.038j", "j1_BB_x_alpha30":[]}
    j1_BB_x_alpha10 = {"M": "1.57 - 0.038j", "j1_BB_x_alpha10":[]}
    
    #gauss bean M = 1.57 - 0.038j
    for xi in xval:
        print("i:", cnt)
        #parou no 37

        j1_PW_x['j1_PW_x'].append({
            "j1z": np.float64(J1(xi, M, epsilonR2(M, 1), 1)),
            "x": xi
        })
        json.dump(j1_PW_x, open(os.path.abspath(path + 'j1_PW_x.json'), 'w'))

        j1_BB_x_alpha1['j1_BB_x_alpha1'].append({
            "j1z": np.float64(I_z(M, xi, epsilonR2(M, 1), 0, math.radians(1), 0, 0)),
            "x": xi
        })
        json.dump(j1_BB_x_alpha1, open(os.path.abspath(path + 'j1_BB_x_alpha1.json'), 'w'))

        j1_BB_x_alpha30['j1_BB_x_alpha30'].append({
            "j1z": np.float64(I_z(M, xi, epsilonR2(M, 1), 0, math.radians(30), 0, 0)),
            "x": xi
        })
        json.dump(j1_BB_x_alpha30, open(os.path.abspath(path + 'j1_BB_x_alpha30.json'), 'w'))

        j1_BB_x_alpha10['j1_BB_x_alpha10'].append({
            "j1z": np.float64(I_z(M, xi, epsilonR2(M, 1), 0, math.radians(10), 0, 0)),
            "x": xi
        })
        json.dump(j1_BB_x_alpha10, open(os.path.abspath(path + 'j1_BB_x_alpha10.json'), 'w'))
        
        cnt = cnt + 1


def cache_BB_Jz_x_5():

    eps = epsilonR2(M, 1)
    eps1 = epsilonR2(M1, 1)
    eps2 = epsilonR2(M2, 1)
    eps3 = epsilonR2(M3, 1)

    cnt = 0
    '''
    jz_BB_x_M = {"jz_BB_x_M_5":[]}
    jz_BB_x_M1 = {"jz_BB_x_M1_5":[]}
    jz_BB_x_M2 = {"jz_BB_x_M2_5":[]}
    jz_BB_x_M3 = {"jz_BB_x_M3_5":[]}
    '''

    with open(os.path.abspath(path + 'jz_BB_x_M_5.json'), 'r') as jz_BB_x_M:
        jz_BB_x_M = json.load(jz_BB_x_M)

    with open(os.path.abspath(path + 'jz_BB_x_M1_5.json'), 'r') as jz_BB_x_M1:
        jz_BB_x_M1 = json.load(jz_BB_x_M1)
    
    with open(os.path.abspath(path + 'jz_BB_x_M2_5.json'), 'r') as jz_BB_x_M2:
        jz_BB_x_M2 = json.load(jz_BB_x_M2)

    with open(os.path.abspath(path + 'jz_BB_x_M3_5.json'), 'r') as jz_BB_x_M3:
        jz_BB_x_M3 = json.load(jz_BB_x_M3)

    xval = np.linspace(0.1, 20, 100)
    #print(xval[42:])

    #calculo do espectro de Iz, ao longo de diferentes tamanhos de particulas, para diferentes indices de refração 
    #parou no 63 -> colocar 64
    xval = xval[64:]
    cnt = 64
    for xi in xval:
        print("i:", cnt)

        jz_BB_x_M['jz_BB_x_M_5'].append({
            "jz": np.float64(I_z(M, xi, eps, 0, math.radians(5), 0, 0)),
            "x": xi
        })
        json.dump(jz_BB_x_M, open(os.path.abspath(path + 'jz_BB_x_M_5.json'), 'w'))

        jz_BB_x_M1['jz_BB_x_M1_5'].append({
            "jz": np.float64(I_z(M1, xi, eps1, 0, math.radians(5), 0, 0)),
            "x": xi
        })
        json.dump(jz_BB_x_M1, open(os.path.abspath(path + 'jz_BB_x_M1_5.json'), 'w'))

        jz_BB_x_M2['jz_BB_x_M2_5'].append({
            "jz": np.float64(I_z(M2, xi, eps2, 0, math.radians(5), 0, 0)),
            "x": xi
        })
        json.dump(jz_BB_x_M2, open(os.path.abspath(path + 'jz_BB_x_M2_5.json'), 'w'))

        jz_BB_x_M3['jz_BB_x_M3_5'].append({
            "jz": np.float64(I_z(M3, xi, eps3, 0, math.radians(5), 0, 0)),
            "x": xi
        })
        json.dump(jz_BB_x_M3, open(os.path.abspath(path + 'jz_BB_x_M3_5.json'), 'w'))
        
        cnt = cnt + 1
        

def cache_BB_Jz_x_30():

    eps = epsilonR2(M, 1)
    eps1 = epsilonR2(M1, 1)
    eps2 = epsilonR2(M2, 1)
    eps3 = epsilonR2(M3, 1)

    cnt = 0

    jz_BB_x_M = {"j1_BB_x_M_30":[]}
    jz_BB_x_M1 = {"j1_BB_x_M1_30":[]}
    jz_BB_x_M2 = {"j1_BB_x_M2_30":[]}
    jz_BB_x_M3 = {"j1_BB_x_M3_30":[]}

    xval = np.linspace(0.1, 20, 100)

    #calculo do espectro de Iz, ao longo de diferentes tamanhos de particulas, para diferentes indices de refração 
    for xi in xval:
        print("i:", cnt)

        jz_BB_x_M['jz_BB_x_M_30'].append({
            "jz": I_z(M, xi, eps, 0, math.radians(30), 0, 0),
            "x": xi
        })
        json.dump(jz_BB_x_M, open(os.path.abspath(path + 'jz_BB_x_M_30.json'), 'w'))

        jz_BB_x_M1['jz_BB_x_M1_30'].append({
            "jz": I_z(M1, xi, eps1, 0, math.radians(30), 0, 0),
            "x": xi
        })
        json.dump(jz_BB_x_M1, open(os.path.abspath(path + 'jz_BB_x_M1_30.json'), 'w'))

        jz_BB_x_M2['jz_BB_x_M2_30'].append({
            "jz": I_z(M2, xi, eps2, 0, math.radians(30), 0, 0),
            "x": xi
        })
        json.dump(jz_BB_x_M2, open(os.path.abspath(path + 'jz_BB_x_M2_30.json'), 'w'))

        jz_BB_x_M3['jz_BB_x_M3_30'].append({
            "jz": I_z(M3, xi, eps3, 0, math.radians(30), 0, 0),
            "x": xi
        })
        json.dump(jz_BB_x_M3, open(os.path.abspath(path + 'jz_BB_x_M3_30.json'), 'w'))
        
        cnt = cnt + 1




def cache_load_BB_Jz_x_5():

    jn1_gauss_blue = []
    jn2_gauss_red = []
    jn3_gauss_green = []
    jn4_gauss_orange = []
    x = []
    x_axis = np.linspace(0, 20, 100)
    count = np.linspace(0, 99, 100, dtype=int)

    with open(os.path.abspath(path + 'jz_BB_x_M_5.json'), 'r') as jz_BB_x_M:
        jz_BB_x_M = json.load(jz_BB_x_M)

    with open(os.path.abspath(path + 'jz_BB_x_M1_5.json'), 'r') as jz_BB_x_M1:
        jz_BB_x_M1 = json.load(jz_BB_x_M1)
    
    with open(os.path.abspath(path + 'jz_BB_x_M2_5.json'), 'r') as jz_BB_x_M2:
        jz_BB_x_M2 = json.load(jz_BB_x_M2)

    with open(os.path.abspath(path + 'jz_BB_x_M3_5.json'), 'r') as jz_BB_x_M3:
        jz_BB_x_M3 = json.load(jz_BB_x_M3)

    for i in count:
        jn1_gauss_blue.append(jz_BB_x_M['jz_BB_x_M_5'][i]['j1'])
        jn2_gauss_red.append(jz_BB_x_M1['jz_BB_x_M1_5'][i]['j1'])
        jn3_gauss_green.append(jz_BB_x_M2['jz_BB_x_M2_5'][i]['j1'])
        jn4_gauss_orange.append(jz_BB_x_M3['jz_BB_x_M3_5'][i]['j1'])
        x.append(jz_BB_x_M3['jz_BB_x_M3_5'][i]["x"])

    plt.figure(figsize=[8,6])
    print(M, M1, M2, M3)
    plt.plot(x, jn4_gauss_orange,"y--", label= r'$M = 1.57 - 1j$' )
    plt.plot(x, jn2_gauss_red,"r", label= r'$M = 1.57- 0.38j$' )
    plt.plot(x, jn1_gauss_blue,"b", label= r'$M = 1.57 - 0.038j$' )
    plt.plot(x, jn3_gauss_green,"g", label= r'$M = 1.57-0.01j$' )
    plt.plot(x, np.linspace(0, 0, 100), "k--")

    plt.xlabel(r' x $(\mu m)$')
    plt.ylabel(r'Assymetry Factor $I_z$')
    plt.title(r'Para $\alpha = 5^{\circ}$')
    plt.legend()
    plt.grid()
    plt.savefig(r'D:/Documentos/USP/TCC/Simulações - TCC - BesselBeam/besselBeam/graphs/jz_BB_x_M3_5.png', format='png', dpi=500)
    plt.show()

def cache_load_BB_Jz_x_30():

    jn1_gauss_blue = []
    jn2_gauss_red = []
    jn3_gauss_green = []
    jn4_gauss_orange = []
    x = []
    x_axis = np.linspace(0, 20, 100)
    count = np.linspace(0, 99, 100, dtype=int)

    with open(os.path.abspath(path + 'jz_BB_x_M_30.json'), 'r') as jz_BB_x_M:
        jz_BB_x_M = json.load(jz_BB_x_M)

    with open(os.path.abspath(path + 'jz_BB_x_M1_30.json'), 'r') as jz_BB_x_M1:
        jz_BB_x_M1 = json.load(jz_BB_x_M1)
    
    with open(os.path.abspath(path + 'jz_BB_x_M2_30.json'), 'r') as jz_BB_x_M2:
        jz_BB_x_M2 = json.load(jz_BB_x_M2)

    with open(os.path.abspath(path + 'jz_BB_x_M3_30.json'), 'r') as jz_BB_x_M3:
        jz_BB_x_M3 = json.load(jz_BB_x_M3)

    for i in count:
        jn1_gauss_blue.append(jz_BB_x_M['jz_BB_x_M_30'][i]['j1'])
        jn2_gauss_red.append(jz_BB_x_M1['jz_BB_x_M1_30'][i]['j1'])
        jn3_gauss_green.append(jz_BB_x_M2['jz_BB_x_M2_30'][i]['j1'])
        jn4_gauss_orange.append(jz_BB_x_M3['jz_BB_x_M3_30'][i]['j1'])
        x.append(jz_BB_x_M3['jz_BB_x_M3_30'][i]["x"])

    plt.figure(figsize=[8,6])
    print(M, M1, M2, M3)
    plt.plot(x, jn4_gauss_orange,"y--", label= r'$M = 1.57 - 1j$' )
    plt.plot(x, jn2_gauss_red,"r", label= r'$M = 1.57- 0.38j$' )
    plt.plot(x, jn1_gauss_blue,"b", label= r'$M = 1.57 - 0.038j$' )
    plt.plot(x, jn3_gauss_green,"g", label= r'$M = 1.57-0.01j$' )
    plt.plot(x, np.linspace(0, 0, 100), "k--")

    plt.xlabel(r' x $(\mu m)$')
    plt.ylabel(r'Assymetry Factor $I_z$')
    plt.title(r'Para $\alpha = 30^{\circ}$')
    plt.legend()
    plt.grid()
    plt.savefig(r'D:/Documentos/USP/TCC/Simulações - TCC - BesselBeam/besselBeam/graphs/BB_Jz_x_30.png', format='png', dpi=500)
    plt.show()



def cache_load_FW1_J1_z():

    jn1_gauss_blue = []
    jn2_gauss_red = []
    jn3_gauss_green = []

    with open(os.path.abspath(path + 'FW_J1_01_z.json'), 'r') as FW_J1_01_z:
        FW_J1_01_z = json.load(FW_J1_01_z)
    with open(os.path.abspath(path + 'FW_J1_3_z.json'), 'r') as FW_J1_3_z:
        FW_J1_3_z = json.load(FW_J1_3_z)
    with open(os.path.abspath(path + 'FW_J1_8_z.json'), 'r') as FW_J1_8_z:
        FW_J1_8_z = json.load(FW_J1_8_z)
    


    z_axis = np.linspace(0, 400, 200)
    count = np.linspace(0, 199, 200, dtype=int)

    #M = 1.57 - 0.038j
    for i in count:

        jn1_gauss_blue.append(FW_J1_01_z['j1_FW_01'][i]['j1']*250)
        jn2_gauss_red.append(FW_J1_3_z['j1_FW_3'][i]['j1'])
        jn3_gauss_green.append(FW_J1_8_z['j1_FW_8'][i]['j1'])


    plt.figure(figsize=[10,8])
    plt.plot(z_axis, jn1_gauss_blue,"b", label= "x = 0.1" )
    plt.plot(z_axis, jn2_gauss_red,"r", label= "x = 3" )    
    plt.plot(z_axis, jn3_gauss_green,"g", label= "x = 8" )

    mpl.rcParams['text.usetex'] = True
    mpl.rcParams['text.latex.preamble'] = r'\usepackage{{amsmath}}'
    plt.xlabel(r'$z_0$ $\mu m$')
    plt.ylabel(r'Assymetry Factor $J_1$')
    plt.title(r"$F(z) = $ $\left\{"
            r"\begin{matrix}"
            r"1, & \text{ se }\frac{3L}{8} \leq z \leq \frac{5L}{8} \\"
            r"0, & c.c"
            r"\end{matrix}"
            r"\right. $")
    plt.legend()
    plt.grid()
    plt.savefig(r'D:/Documentos/USP/TCC/Simulações - TCC - BesselBeam/besselBeam/graphs/FW1_J1_z.png', format='png', dpi=500)
    plt.show()


def cache_load_FW2_J1_z():

    jn1_gauss_blue = []
    jn2_gauss_red = []
    jn3_gauss_green = []

    with open(os.path.abspath(path + 'FW2_J1_01_z.json'), 'r') as FW2_J1_01_z:
        FW2_J1_01_z = json.load(FW2_J1_01_z)
    with open(os.path.abspath(path + 'FW2_J1_3_z.json'), 'r') as FW2_J1_3_z:
        FW2_J1_3_z = json.load(FW2_J1_3_z)
    with open(os.path.abspath(path + 'FW2_J1_8_z.json'), 'r') as FW2_J1_8_z:
        FW2_J1_8_z = json.load(FW2_J1_8_z)
    

    z_axis = np.linspace(0, 400, 200)
    count = np.linspace(0, 199, 200, dtype=int)

    #M = 1.57 - 0.038j
    for i in count:

        jn1_gauss_blue.append(FW2_J1_01_z['j1_FW2_01'][i]['j1']*250)
        jn2_gauss_red.append(FW2_J1_3_z['j1_FW2_3'][i]['j1'])
        jn3_gauss_green.append(FW2_J1_8_z['j1_FW2_8'][i]['j1'])


    plt.figure(figsize=[10,8])
    plt.plot(z_axis, jn1_gauss_blue,"b", label= "x = 0.1" )
    plt.plot(z_axis, jn2_gauss_red,"r", label= "x = 3" )    
    plt.plot(z_axis, jn3_gauss_green,"g", label= "x = 8" )

    mpl.rcParams['text.usetex'] = True
    mpl.rcParams['text.latex.preamble'] = r'\usepackage{{amsmath}}'
    plt.xlabel(r'$z_0$ $\mu m$')
    plt.ylabel(r'Assymetry Factor $J_1$')
    plt.title(r"$F(z) = $ $\left\{"
            r"\begin{matrix}"
            r"e^{10^4 (z - 6L/8)/2}, & \text{ se }\frac{2L}{8} \leq z \leq \frac{6L}{8} \\"
            r"0, & \text{ c.c }"
            r"\end{matrix}"
            r"\right. $")
    plt.legend()
    plt.grid()
    plt.savefig(r'D:/Documentos/USP/TCC/Simulações - TCC - BesselBeam/besselBeam/graphs/FW2_J1_z.png', format='png', dpi=500)
    plt.show()


def cache_load_FW3_J1_z():

    jn1_gauss_blue = []
    jn2_gauss_red = []
    jn3_gauss_green = []

    with open(os.path.abspath(path + 'FW3_J1_01_z.json'), 'r') as FW3_J1_01_z:
        FW3_J1_01_z = json.load(FW3_J1_01_z)
    with open(os.path.abspath(path + 'FW3_J1_3_z.json'), 'r') as FW3_J1_3_z:
        FW3_J1_3_z = json.load(FW3_J1_3_z)
    with open(os.path.abspath(path + 'FW3_J1_8_z.json'), 'r') as FW3_J1_8_z:
        FW3_J1_8_z = json.load(FW3_J1_8_z)
    

    z_axis = np.linspace(0, 400, 200)
    count = np.linspace(0, 199, 200, dtype=int)

    #M = 1.57 - 0.038j
    for i in count:

        jn1_gauss_blue.append(FW3_J1_01_z['j1_FW3_01'][i]['j1']*250)
        jn2_gauss_red.append(FW3_J1_3_z['j1_FW3_3'][i]['j1'])
        jn3_gauss_green.append(FW3_J1_8_z['j1_FW3_8'][i]['j1'])


    plt.figure(figsize=[10,8])
    plt.plot(z_axis, jn1_gauss_blue,"b", label= "x = 0.1" )
    plt.plot(z_axis, jn2_gauss_red,"r", label= "x = 3" )    
    plt.plot(z_axis, jn3_gauss_green,"g", label= "x = 8" )

    mpl.rcParams['text.usetex'] = True
    mpl.rcParams['text.latex.preamble'] = r'\usepackage{{amsmath}}'
    plt.xlabel(r'$z_0$ $\mu m$')
    plt.ylabel(r'Assymetry Factor $J_1$')
    plt.title(r"$F(z) = $ $\left\{"
            r"\begin{matrix}"
            r"exp \left[ -5 \frac{(z - 0.5 L)^2}{L^2}\right] cos\left(\frac{6 \pi z}{L}\right), & \text{ se }\frac{L}{12} \leq z \leq \frac{11L}{12} \\"
            r"0, & c.c"
            r"\end{matrix}"
            r"\right. $")
    plt.legend()
    plt.grid()
    plt.savefig(r'D:/Documentos/USP/TCC/Simulações - TCC - BesselBeam/besselBeam/graphs/FW3_J1_z.png', format='png', dpi=500)
    plt.show()

def cache_load_FW4_J1_z():

    jn1_gauss_blue = []
    jn2_gauss_red = []
    jn3_gauss_green = []

    with open(os.path.abspath(path + 'FW4_J1_01_z.json'), 'r') as FW4_J1_01_z:
        FW4_J1_01_z = json.load(FW4_J1_01_z)
    with open(os.path.abspath(path + 'FW4_J1_3_z.json'), 'r') as FW4_J1_3_z:
        FW4_J1_3_z = json.load(FW4_J1_3_z)
    with open(os.path.abspath(path + 'FW4_J1_8_z.json'), 'r') as FW4_J1_8_z:
        FW4_J1_8_z = json.load(FW4_J1_8_z)
    

    z_axis = np.linspace(0, 400, 200)
    count = np.linspace(0, 199, 200, dtype=int)

    #M = 1.57 - 0.038j
    for i in count:

        jn1_gauss_blue.append(FW4_J1_01_z['j1_FW4_01'][i]['j1']*250)
        jn2_gauss_red.append(FW4_J1_3_z['j1_FW4_3'][i]['j1'])
        jn3_gauss_green.append(FW4_J1_8_z['j1_FW4_8'][i]['j1'])


    plt.figure(figsize=[10,8])
    plt.plot(z_axis, jn1_gauss_blue,"b", label= "x = 0.1" )
    plt.plot(z_axis, jn2_gauss_red,"r", label= "x = 3" )    
    plt.plot(z_axis, jn3_gauss_green,"g", label= "x = 8" )

    mpl.rcParams['text.usetex'] = True
    mpl.rcParams['text.latex.preamble'] = r'\usepackage{{amsmath}}'
    plt.xlabel(r'$z_0$ $\mu m$')
    plt.ylabel(r'Assymetry Factor $J_1$')

    plt.legend()
    plt.grid()
    plt.savefig(r'D:/Documentos/USP/TCC/Simulações - TCC - BesselBeam/besselBeam/graphs/FW4_J1_z.png', format='png', dpi=500)
    plt.show()



def cache_load_FW1_J1_x():

    jn1_gauss_blue = []
    jn2_gauss_red = []
    x = []
    x_axis = np.linspace(0, 20, 100)
    count = np.linspace(0, 99, 100, dtype=int)

    with open(os.path.abspath(path + 'FW1_j1_x_L2.json'), 'r') as j1_FW1_x_L2:
        j1_FW1_x_L2 = json.load(j1_FW1_x_L2)

    with open(os.path.abspath(path + 'FW1_j1_x_L4.json'), 'r') as j1_FW1_x_L4:
        j1_FW1_x_L4 = json.load(j1_FW1_x_L4)

    #M = 1.57 - 0.038j
    for i in count:
        jn1_gauss_blue.append(j1_FW1_x_L2['j1_FW1_x_L2'][i]['j1'])
        jn2_gauss_red.append(j1_FW1_x_L4['j1_FW1_x_L4'][i]['j1'] * 25000)
        x.append(j1_FW1_x_L2['j1_FW1_x_L2'][i]["x"])

    plt.figure(figsize=[10,8])
    plt.plot(x, jn1_gauss_blue,"b", label= r'$z = \frac{L}{2}$')
    plt.plot(x, jn2_gauss_red,"r", label= r'$z = \frac{L}{4} \times 25000$' )    
    plt.xlabel(r'$x$ $\mu m$')
    plt.ylabel(r'Assymetry Factor $J_1$')
    plt.legend()
    plt.grid()
    plt.savefig(r'D:/Documentos/USP/TCC/Simulações - TCC - BesselBeam/besselBeam/graphs/FW1_J1_x.png', format='png', dpi=500)
    plt.show()

def cache_load_FW3_J1_x():

    jn1_gauss_blue = []
    jn2_gauss_red = []
    x = []
    x_axis = np.linspace(0, 20, 100)
    count = np.linspace(0, 99, 100, dtype=int)

    with open(os.path.abspath(path + 'FW3_j1_x_L2.json'), 'r') as j1_FW3_x_L2:
        j1_FW3_x_L2 = json.load(j1_FW3_x_L2)

    with open(os.path.abspath(path + 'FW3_j1_x_L4.json'), 'r') as j1_FW3_x_L4:
        j1_FW3_x_L4 = json.load(j1_FW3_x_L4)

    #M = 1.57 - 0.038j
    for i in count:
        jn1_gauss_blue.append(j1_FW3_x_L2['j1_FW3_x_L2'][i]['j1'])
        jn2_gauss_red.append(j1_FW3_x_L4['j1_FW3_x_L4'][i]['j1'] * 100)
        x.append(j1_FW3_x_L4['j1_FW3_x_L4'][i]["x"])


    plt.figure(figsize=[10,8])
    plt.plot(x, jn1_gauss_blue,"b", label= r'$z = \frac{L}{2}$')
    plt.plot(x, jn2_gauss_red,"r", label= r'$z = \frac{L}{4} \times 100$' )    
    plt.xlabel(r'$x$ $\mu m$')
    plt.ylabel(r'Assymetry Factor $J_1$')
    plt.legend()
    plt.grid()
    plt.savefig(r'D:/Documentos/USP/TCC/Simulações - TCC - BesselBeam/besselBeam/graphs/FW3_J1_x.png', format='png', dpi=500)
    plt.show()

def cache_load_FW2_J1_x():

    jn1_gauss_blue = []
    jn2_gauss_red = []
    x = []
    x_axis = np.linspace(0, 20, 100)
    count = np.linspace(0, 99, 100, dtype=int)

    with open(os.path.abspath(path + 'j1_FW2_x_L2.json'), 'r') as FW2_j1_x_L2:
        FW2_j1_x_L2 = json.load(FW2_j1_x_L2)

    with open(os.path.abspath(path + 'j1_FW2_x_L4.json'), 'r') as FW2_j1_x_L4:
        FW2_j1_x_L4 = json.load(FW2_j1_x_L4)

    #M = 1.57 - 0.038j
    for i in count:
        jn1_gauss_blue.append(FW2_j1_x_L2['j1_FW2_x_L2'][i]['j1'])
        jn2_gauss_red.append(FW2_j1_x_L4['j1_FW2_x_L4'][i]['j1'] )
        x.append(FW2_j1_x_L4['j1_FW2_x_L4'][i]["x"])


    plt.figure(figsize=[10,8])
    plt.plot(x, jn1_gauss_blue,"b", label= r'$z = \frac{L}{2}$')
    plt.plot(x, jn2_gauss_red,"r", label= r'$z = \frac{L}{4} $' )    
    plt.xlabel(r'$x$ $\mu m$')
    plt.ylabel(r'Assymetry Factor $J_1$')
    plt.legend()
    plt.grid()
    plt.savefig(r'D:/Documentos/USP/TCC/Simulações - TCC - BesselBeam/besselBeam/graphs/FW2_J1_x.png', format='png', dpi=500)
    plt.show()

def cache_load_FW1_Psi():

    Psi = []
    with open(os.path.abspath(path + 'FW1_Psi.json'), 'r') as FW1_Psi:
        FW1_Psi = json.load(FW1_Psi)


    z_axis = np.linspace(0, 400, 200)
    count = np.linspace(0, 199, 200, dtype=int)

    #M = 1.57 - 0.038j
    for i in count:
        Psi.append(FW1_Psi['FW1_Psi'][i]['Psi'])


    plt.figure(figsize=[10,8])
    plt.plot(z_axis, Psi,"b")

    plt.title(r'$|\Psi (0, z)|^2$')
    plt.xlabel(r'$z (\mu m)$')
    plt.grid()
    plt.savefig(r'D:/Documentos/USP/TCC/Simulações - TCC - BesselBeam/besselBeam/graphs/FW1_Psi.png', format='png', dpi=500)
    plt.show()

def cache_load_FW2_Psi():

    Psi = []
    with open(os.path.abspath(path + 'FW2_Psi.json'), 'r') as FW2_Psi:
        FW2_Psi = json.load(FW2_Psi)


    z_axis = np.linspace(0, 400, 200)
    count = np.linspace(0, 199, 200, dtype=int)

    #M = 1.57 - 0.038j
    for i in count:
        Psi.append(FW2_Psi['FW2_Psi'][i]['Psi'])


    plt.figure(figsize=[10,8])
    plt.plot(z_axis, Psi,"b")

    plt.title(r'$|\Psi (0, z)|^2$')
    plt.xlabel(r'$z (\mu m)$')
    plt.grid()
    plt.savefig(r'D:/Documentos/USP/TCC/Simulações - TCC - BesselBeam/besselBeam/graphs/FW2_Psi.png', format='png', dpi=500)
    plt.show()

def cache_load_FW3_Psi():

    Psi = []
    with open(os.path.abspath(path + 'FW3_Psi.json'), 'r') as FW3_Psi:
        FW3_Psi = json.load(FW3_Psi)


    z_axis = np.linspace(0, 400, 200)
    count = np.linspace(0, 199, 200, dtype=int)

    #M = 1.57 - 0.038j
    for i in count:
        Psi.append(FW3_Psi['FW3_Psi'][i]['Psi'])


    plt.figure(figsize=[10,8])
    plt.plot(z_axis, Psi,"b")

    plt.title(r'$|\Psi (0, z)|^2$')
    plt.xlabel(r'$z (\mu m)$')
    plt.grid()
    plt.savefig(r'D:/Documentos/USP/TCC/Simulações - TCC - BesselBeam/besselBeam/graphs/FW3_Psi.png', format='png', dpi=500)
    plt.show()

def cache_load_FW4_Psi():

    Psi = []
    with open(os.path.abspath(path + 'FW4_Psi.json'), 'r') as FW3_Psi:
        FW3_Psi = json.load(FW3_Psi)


    z_axis = np.linspace(0, 400, 200)
    count = np.linspace(0, 199, 200, dtype=int)

    #M = 1.57 - 0.038j
    for i in count:
        Psi.append(FW3_Psi['FW4_Psi'][i]['Psi'])


    plt.figure(figsize=[10,8])
    plt.plot(z_axis, Psi,"b")

    plt.title(r'$|\Psi (0, z)|^2$')
    plt.xlabel(r'$z (\mu m)$')
    plt.grid()
    plt.savefig(r'D:/Documentos/USP/TCC/Simulações - TCC - BesselBeam/besselBeam/graphs/FW4_Psi.png', format='png', dpi=500)
    plt.show()


def cache_load_BB_J1_x_alphas():

    jn1_gauss_blue = []
    jn2_gauss_red = []
    jn3_gauss_green = []
    jn4_gauss_orange = []
    x = []
    x_axis = np.linspace(0, 20, 100)
    count = np.linspace(0, 99, 100, dtype=int)

    with open(os.path.abspath(path + 'j1_PW_x.json'), 'r') as j1_PW_x:
        j1_PW_x = json.load(j1_PW_x)

    with open(os.path.abspath(path + 'j1_BB_x_alpha1.json'), 'r') as j1_BB_x_alpha1:
        j1_BB_x_alpha1 = json.load(j1_BB_x_alpha1)

    with open(os.path.abspath(path + 'j1_BB_x_alpha10.json'), 'r') as j1_BB_x_alpha10:
        j1_BB_x_alpha10 = json.load(j1_BB_x_alpha10)
    
    with open(os.path.abspath(path + 'j1_BB_x_alpha30.json'), 'r') as j1_BB_x_alpha30:
        j1_BB_x_alpha30 = json.load(j1_BB_x_alpha30)

    '''
    with open(os.path.abspath(path + 'j1_BB_x_alpha60.json'), 'r') as j1_BB_x_alpha60:
        j1_BB_x_alpha60 = json.load(j1_BB_x_alpha60)
    '''

    for i in count:
        jn1_gauss_blue.append(j1_PW_x['j1_PW_x'][i]['j1z'])
        jn2_gauss_red.append(j1_BB_x_alpha1['j1_BB_x_alpha1'][i]['j1z'])
        jn3_gauss_green.append(j1_BB_x_alpha30['j1_BB_x_alpha30'][i]['j1z'])
        jn4_gauss_orange.append(j1_BB_x_alpha10['j1_BB_x_alpha10'][i]['j1z'])
        x.append(j1_BB_x_alpha10['j1_BB_x_alpha60'][i]["x"])

    plt.figure(figsize=[8,6])
    print(M, M1, M2, M3)
    plt.plot(x, jn1_gauss_blue,"k", label= r'Plane Wave' )
    plt.plot(x, jn2_gauss_red,"r-.", label= r'$\alpha = 1^{\circ}$' )
    plt.plot(x, jn3_gauss_green,"b--", label= r'$\alpha = 30^{\circ}$' )
    plt.plot(x, jn4_gauss_orange,"g", label= r'$\alpha = 10^{\circ}$' )
    #plt.plot(x, np.linspace(0, 0, 100), "k--")

    plt.xlabel(r'Size Paramter x $(\mu m)$')
    plt.ylabel(r'Assymetry Factor $J_{1, z}$')
    #plt.title(r'Para $\alpha = 5^{\circ}$')
    plt.legend()
    plt.grid()
    plt.savefig(r'D:/Documentos/USP/TCC/Simulações - TCC - BesselBeam/besselBeam/graphs/BB_J1_x_alphas.png', format='png', dpi=500)
    plt.show()

#rodar pontos 
#cache_FW1_J1_z()
#cache_FW2_J1_z()
#cache_FW3_J1_z()

#  salvar pontos - intensidade de campo em z
#cache_FW1_Psi()
#cache_FW2_Psi()
#cache_FW3_Psi()
#cache_FW4_Psi()

# salvar pontos - J1 por x 
#cache_FW1_J1_x()
#cache_FW3_J1_x()
#cache_FW2_J1_x()
#cache_BB_J1_x_alphas()
#cache_BB_Jz_x_5()
#cache_BB_Jz_x_30() #rodar amanha 

#carregar pontos 
#cache_load_FW1_J1_z()
#cache_load_FW2_J1_z()
#cache_load_FW3_J1_z()


#loads de intensidade de campo 
#cache_load_FW1_Psi()
#cache_load_FW2_Psi()
#cache_load_FW3_Psi()
#cache_load_FW4_Psi()

#loads - J1 por x   
#cache_load_FW1_J1_x()
#cache_load_FW2_J1_x()
#cache_load_FW3_J1_x() 
#cache_load_BB_Jz_x_5()
#cache_load_BB_Jz_x_30() 
#cache_load_BB_J1_x_alphas()