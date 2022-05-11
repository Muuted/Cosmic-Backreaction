import scipy
import scipy.integrate

import matplotlib.pyplot as plt

from scipy.integrate import quad

from scipy.integrate import solve_ivp

import numpy as np

'''
r = List[0]
    RR = List[1]
    EE = List[2]
    MM = List[3]
    dMMdr = List[4]
    dRRdr = List[5]
    dEEdr = List[6]
    G = List[7]
    rho_FLRW = List[8]
    r_b = List[9]
    n = List[10]
    m = List[11]
    A = List[12]
    H = List[13]
    Lamb = List[14]
    c = List[15]
'''

def func_E(r, RR, EE, MM, dMMdr, dRRdr, dEEdr, G, rho_c0, r_b, n, m, A, H, Lamb,c):
    """
    r = List[0]
    r_b = List[9]
    n = List[10]
    m = List[11]
    A = List[12]
    """
    if r <= r_b:

        EE = A * r ** 2 * ((r / r_b) ** n - 1) ** m

    else:

        EE = 0

    return EE


def func_dEdr(r, RR, EE, MM, dMMdr, dRRdr, dEEdr, G, rho_c0, r_b, n, m, A, H, Lamb,c):
    """r = List[0]
    r_b = List[9]
    n = List[10]
    m = List[11]
    A = List[12]"""
    
    dEEdr1 = 2 * A * r * ((r / r_b) ** n - 1) ** m

    dEEdr2 = m * n * A * r ** (n + 1) * ((r / r_b) ** n - 1) ** (m - 1)

    dEEdr = dEEdr1 + dEEdr2
    return dEEdr


def func_rho(r, RR, EE, MM, dMMdr, dRRdr, dEEdr, G, rho_c0, r_b, n, m, A, H, Lamb,c):
    """RR = List[1]
    dMMdr = List[4]
    dRRdr = List[5]
    G = List[7]"""
    
    kappa = (8*np.pi*G)/(c**4)*rho_c0

    rho = (2/kappa)*(dMMdr/(RR**2*dRRdr))

    return rho

def func_dMMdr(r, RR, EE, MM, dMMdr, dRRdr, dEEdr, G, rho_c0, r_b, n, m, A, H, Lamb,c):

    const = (4*np.pi*G)/(c**4)

    dMMdr = const * RR**2 * dRRdr

    return dMMdr


def func_dRdt(r, RR, EE, MM, dMMdr, dRRdr, dEEdr, G, rho_FLRW, r_b, n, m, A, H, Lamb,c):
    '''
    RR = List[1]
    EE = List[2]
    MM = List[3]
    Lamb = List[14]'''

    dRRdt = np.sqrt(2 * MM / RR + 2 * EE + (Lamb / 3) * RR ** 2)

    return dRRdt


def func_dRdrdt(r, RR, EE, MM, dMMdr, dRRdr, dEEdr, G, rho_FLRW, r_b, n, m, A, H, Lamb,c):
    
    '''
    RR = List[1]
    EE = List[2]
    MM = List[3]
    dMMdr = List[4]
    dRRdr = List[5]
    dEEdr = List[6]
    Lamb = List[14]
    '''
    
    sqrts = np.sqrt(    2*MM/RR +2*EE + (Lamb/3)*RR**2  )

    extra = 2*dMMdr/RR - 4*MM*dRRdr/RR**2 +2*dEEdr + 2*Lamb*RR*dRRdr/3

    return extra/sqrts


def func_LTB_dSdt(S,t,p):
    RR, dRRdr, MM, rho = S
    #r, EE,dEEdr, dMMdr, dEEdr, G, rho_c0, r_b, n, m, A, H, Lamb,c = p
    r, EE, dEEdr, dMMdr, G, rho_c0, r_b, n, m, A, H_0, Lamb,c = p

    dRdt = func_dRdt(r, RR, EE, MM, dMMdr, dRRdr, dEEdr, G, rho_c0, r_b, n, m, A, H_0, Lamb, c)
    dRdrdt =  func_dRdrdt(r, RR, EE, MM, dMMdr, dRRdr, dEEdr, G, rho_c0, r_b, n, m, A, H_0, Lamb,c)
    dMMdr = func_dMMdr(r, RR, EE, MM, dMMdr, dRRdr, dEEdr, G, rho_c0, r_b, n, m, A, H_0, Lamb, c)
    rho = func_rho(r, RR, EE, MM, dMMdr, dRRdr, dEEdr, G, rho_c0, r_b, n, m, A, H_0, Lamb, c)

    return_list = [ dRdt, dRdrdt, dMMdr, rho ]
    return  return_list


"""
def M(List):
    
    RR = List[1]
    EE = List[2]

    rho_FLRW = List[8]
    H = List[13]
    
    MM =( 4 * np.pi * RR ** 3 * rho_FLRW / 3) * (1 - 2 * EE / (5 * H ** 2 * RR ** 2))

    return MM

def dMdr(List):
    
    RR = List[1]
    EE = List[2]
    dRRdr = List[5]
    dEEdr = List[6]
    rho_FLRW = List[8]
    H = List[13]

    dMMdr1 = 12 * np.pi * RR ** 2 * dRRdr * rho_FLRW / 3

    dMMdr2 = -(8 * np.pi * rho_FLRW / 15 * H ** 2) * (EE * dRRdr + RR * dEEdr)

    return dMMdr1 + dMMdr2

def func_rho(List):
    
    RR = List[1]
    dMMdr = List[4]
    dRRdr = List[5]
    G = List[7]

    rrho = (1 / (4 * np.pi * G)) * dMMdr / (RR ** 2 * dRRdr)


    return rrho  

def dSdt_dRdrdt(S,t,p):
    R, dRdr = S
    r, RR, EE, MM, dMMdr, dRRdr, dEEdr, G, rho_FLRW, r_b, n, m, A, H, Lamb = p
    
    ans_dRdr = ddRdrdt(r, RR, EE, MM, dMMdr, dRRdr, dEEdr, G, rho_FLRW, r_b, n, m, A, H, Lamb)
    ans_R = dRdt(r, RR, EE, MM, dMMdr, dRRdr, dEEdr, G, rho_FLRW, r_b, n, m, A, H, Lamb)
    
    return [ans_R, ans_dRdr]

def dSdt_dRdrdt_ivp(t,S,r, RR, EE, MM, dMMdr, dRRdr, dEEdr, G, rho_FLRW, r_b, n, m, A, H, Lamb):
    R, dRdr = S
    #r, RR, EE, MM, dMMdr, dRRdr, dEEdr, G, rho_FLRW, r_b, n, m, A, H, Lamb = p
    
    ans_dRdr = ddRdrdt(r, RR, EE, MM, dMMdr, dRRdr, dEEdr, G, rho_FLRW, r_b, n, m, A, H, Lamb)
    ans_R = dRdt(r, RR, EE, MM, dMMdr, dRRdr, dEEdr, G, rho_FLRW, r_b, n, m, A, H, Lamb)
    
    return [ans_R, ans_dRdr]


 """
# For the de Sitter model universe, for comparison with, our model of the universe. 

def func_FLRW_R_dRdr(t_vec,r_func,H_0):

    # R = a(t)*r
    # a(t) = (t/t_0)^(3/2)
    t_0 = 2/(3*H_0)

    List = [[],[]]

    for i in range(0,len(t_vec)):

        a_FLRW = (t_vec[i]/t_0)**(3/2)
        R_FLRW = a_FLRW*r_func

        List[0].append(R_FLRW)
        List[1].append(a_FLRW)

    return List
