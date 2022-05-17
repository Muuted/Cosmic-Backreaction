import scipy
import scipy.integrate

import matplotlib.pyplot as plt

from scipy.integrate import quad

from scipy.integrate import solve_ivp

import numpy as np
from The_constants import *


Lamb, A, r_b, n, m, H_0, H_i, G, rho_c0, rho_i0, a_i, t_i, t_0,c= func_constants()

def func_E(r, RR, EE, MM, dMMdr, dRRdr, dEEdr, G, rho_c0, r_b, n, m, A, H, Lamb,c):

    if r <= r_b:

        EE = A * r ** 2 * ((r / r_b) ** n - 1) ** m

    else:

        EE = 0

    return EE


def func_dEdr(r, RR, EE, MM, dMMdr, dRRdr, dEEdr, G, rho_c0, r_b, n, m, A, H, Lamb,c):

    if r <= r_b :
        dEEdr1 = 2 * A * r * ((r / r_b) ** n - 1) ** m

        dEEdr2 = m * n * A * r ** (n + 1) * ((r / r_b) ** n - 1) ** (m - 1)

        dEEdr = dEEdr1 + dEEdr2
    else:
        dEEdr = 0
    return dEEdr


def func_rho(r, RR, EE, MM, dMMdr, dRRdr, dEEdr, G, rho_c0, r_b, n, m, A, H, Lamb,c):

    kappa = (8*np.pi)*rho_c0

    rho = (2/kappa)*(dMMdr/(RR**2*dRRdr))

    return rho

def func_M(r, RR, EE, MM, dMMdr, dRRdr, dEEdr, G, rho_c0, r_b, n, m, A, H, Lamb,c):
    
    M_bg = (4*np.pi*G)/(3*c**2)*r**3*rho_c0
    ana_M = M_bg*1+ dMMdr

    return ana_M

def func_dMMdr(r, RR, EE, MM, dMMdr, dRRdr, dEEdr, G, rho_c0, r_b, n, m, A, H, Lamb,c):
    Lamb, A, r_b, n, m, H_0, H_i, G, rho_c0, rho_i0, a_i, t_i, t_0,c= func_constants()
    """
    our origian function 
    const = (4*np.pi*G)/(c**4)*rho_c0
    dMMdr = const * RR**2 * dRRdr 
    """
    # new function
    #M_bg = (4*np.pi*G)/(3*c**2)*r**3*rho_c0
    #dMMdr = M_bg*((6*EE)/(5*r**2))*(c/(a_i*H_i))**2

    dMMdr = 4 * np.pi * rho_c0*r**2  - (24*np.pi*c**2)/(15*a_i**2*H_i**2)*(EE + r*dEEdr)
    return dMMdr


def func_dRdt(r, RR, EE, MM, dMMdr, dRRdr, dEEdr, G, rho_FLRW, r_b, n, m, A, H, Lamb,c):

    dRRdt = c*np.sqrt(
        2 * G*MM /(c**2* RR) + 2 * EE #+ (Lamb / 3) * RR ** 2
        )

    return dRRdt


def func_dRdrdt(r, RR, EE, MM, dMMdr, dRRdr, dEEdr, G, rho_FLRW, r_b, n, m, A, H, Lamb,c):
    #c = 1
    sqrts = 2*np.sqrt(
        2*G*MM/(c**2*RR) +2*EE #+ (Lamb/3)*RR**2 
        )

    extra = (G)/(c)*(dMMdr/RR - MM*dRRdr/(RR**2)) + c*dEEdr #+ Lamb*RR*dRRdr/3

    return extra/sqrts


def func_LTB_dSdt(S,t,p):
    RR, dRRdr, rho = S
    
    r, EE, dEEdr, MM,dMMdr, G, rho_c0, r_b, n, m, A, H_0, Lamb,c = p

    dRdt = func_dRdt(r, RR, EE, MM, dMMdr, dRRdr, dEEdr, G, rho_c0, r_b, n, m, A, H_0, Lamb, c)
    dRdrdt =  func_dRdrdt(r, RR, EE, MM, dMMdr, dRRdr, dEEdr, G, rho_c0, r_b, n, m, A, H_0, Lamb,c)
    #dMMdr = func_dMMdr(r, RR, EE, MM, dMMdr, dRRdr, dEEdr, G, rho_c0 , r_b, n, m, A, H_0, Lamb, c)
    rho = func_rho(r, RR, EE, MM, dMMdr, dRRdr, dEEdr, G, rho_c0, r_b, n, m, A, H_0, Lamb, c)
    #ana_MM = func_M(r, RR, EE, MM, dMMdr, dRRdr, dEEdr, G, rho_c0, r_b, n, m, A, H_0, Lamb, c)

    return_list = [ dRdt, dRdrdt, rho ]#, ana_MM]
    return  return_list



H0 = 0.07158985 # Hubble sonstant 70 km/s/Mpc in unit 1/Gyr
t0 = 2.0/3/H0 # assume EdS background
c = 306.60139 # speed of light in vacuum in units Mpc/Gyr
M_sun = 1.98855e30
Gyr = 3.15576e16
Mpc = 3.0856778e22
r_b = 40
k_max = 1.3e-7
a_i = 1.0/1100
G = 6.673e-11*M_sun*Gyr*Gyr/Mpc/Mpc/Mpc

rho_ic = 3.0*H0*H0/(8.0*np.pi*G)/a_i**3
H_i = np.sqrt( 8*np.pi*G*rho_ic/3.0 )

#print(rho_ic)
#print(H0/H_i)

