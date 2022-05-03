import scipy.integrate

import matplotlib.pyplot as plt

from scipy.integrate import quad

from scipy.integrate import solve_ivp

import numpy as np



def E(r, RR, EE, MM, dMMdr, dRRdr, dEEdr, G, rho_FLRW, r_b, n, m, A, H, Lamb):

    if r <= r_b:

        EE = A * r ** 2 * ((r / r_b) ** n - 1) ** m

    else:

        EE = 0

    return EE



def M(r, RR, EE, MM, dMMdr, dRRdr, dEEdr, G, rho_FLRW, r_b, n, m, A, H, Lamb):


    MM =( 4 * np.pi * RR ** 3 * rho_FLRW / 3) * (1 - 2 * EE / (5 * H ** 2 * RR ** 2))

    return MM



def dEdr(r, RR, EE, MM, dMMdr, dRRdr, dEEdr, G, rho_FLRW, r_b, n, m, A, H, Lamb):

    dEEdr1 = 2 * A * r * ((r / r_b) ** n - 1) ** m

    dEEdr2 = m * n * A * r ** (n + 1) * ((r / r_b) ** n - 1) ** (m - 1)

    return dEEdr1 + dEEdr2



def dMdr(r, RR, EE, MM, dMMdr, dRRdr, dEEdr, G, rho_FLRW, r_b, n, m, A, H, Lamb):

    dMMdr1 = 12 * np.pi * RR ** 2 * dRRdr * rho_FLRW / 3

    dMMdr2 = -(8 * np.pi * rho_FLRW / 15 * H ** 2) * (EE * dRRdr + RR * dEEdr)

    return dMMdr1 + dMMdr2



def dRdt(r, RR, EE, MM, dMMdr, dRRdr, dEEdr, G, rho_FLRW, r_b, n, m, A, H, Lamb):

    dRRdt = np.sqrt(2 * MM / RR + 2 * EE + (Lamb / 3) * RR ** 2)

    return dRRdt



def rho(r, RR, EE, MM, dMMdr, dRRdr, dEEdr, G, rho_FLRW, r_b, n, m, A, H, Lamb):

    rrho = (1 / (4 * np.pi * G)) * dMMdr / (RR ** 2 * dRRdr)

    return rrho


def ddRdrdt(r, RR, EE, MM, dMMdr, dRRdr, dEEdr, G, rho_FLRW, r_b, n, m, A, H, Lamb):

    sqrts = np.sqrt(

        2*MM/RR +2*EE + (Lamb/3)*RR**2
    )


    extra = 2*dMMdr/RR - 4*MM*dRRdr/RR**2 +2*dEEdr + 2*Lamb*RR*dRRdr/3


    return extra/sqrts


'''
def units_converter(number,start_units,end_unitsnu

    one_Mpc = 3.086e22 # m
    one_Gy = 3.156e16 #s
    one_solar_mass = 1.989e30 #kg
    


    one_Mpc = 3.086e22 # m
    one_Gy = 3.156e16 #s
    one_solar_mass = 1.989e30 #kg
    '''