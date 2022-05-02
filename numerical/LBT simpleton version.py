from re import S
import scipy.integrate
import matplotlib.pyplot as plt
from scipy.integrate import quad
from scipy.integrate import solve_ivp
import numpy as np


def E(r, RR, EE, MM, dMMdr, dRRdr, dEEdr, G, rho_i, r_b, n, m, A, H, Lamb):
    if r <= r_b:
        EE = A * r ** 2 * ((r / r_b) ** n - 1) ** m
    else:
        EE = 0
    return EE


def M(r, RR, EE, MM, dMMdr, dRRdr, dEEdr, G, rho_i, r_b, n, m, A, H, Lamb):
    MM = 4 * np.pi * RR ** 3 * rho_i / 3 * (1 - 2 * EE / (5 * H ** 2 * RR ** 2))
    return MM


def dEdr(r, RR, EE, MM, dMMdr, dRRdr, dEEdr, G, rho_i, r_b, n, m, A, H, Lamb):
    dEEdr1 = 2 * A * r * ((r / r_b) ** n - 1) ** m
    dEEdr2 = m * n * A * r ** (n + 1) * ((r / r_b) ** n - 1) ** (m - 1)
    return dEEdr1 + dEEdr2


def dMdr(r, RR, EE, MM, dMMdr, dRRdr, dEEdr, G, rho_i, r_b, n, m, A, H, Lamb):
    dMMdr1 = 12 * np.pi * RR ** 2 * dRRdr * rho_i / 3
    dMMdr2 = -(8 * np.pi * rho_i / 15 * H ** 2) * (EE * dRRdr + RR * dEEdr)
    return dMMdr1 + dMMdr2


def dRdt(r, RR, EE, MM, dMMdr, dRRdr, dEEdr, G, rho_i, r_b, n, m, A, H, Lamb):
    dRRdt = np.sqrt(2 * MM / RR + 2 * EE + (Lamb / 3) * RR ** 2)
    return dRRdt


def rho(r, RR, EE, MM, dMMdr, dRRdr, dEEdr, G, rho_i, r_b, n, m, A, H, Lamb):
    rrho = (1 / (4 * np.pi * G)) * dMMdr / (RR ** 2 * dRRdr)
    return rrho

def ddRdrdt(r, RR, EE, MM, dMMdr, dRRdr, dEEdr, G, rho_i, r_b, n, m, A, H, Lamb):
    sqrts = np.sqrt(
        2*MM/RR +2*EE + (Lamb/3)*RR**2
    )

    extra = 2*dMMdr/RR - 4*MM*dRRdr/RR**2 +2*dEEdr + 2*Lamb*RR*dRRdr/3

    return extra/sqrts


'''
Making a for loop on the radius, where we will integrate with the ODEint; on time;
for each loop
'''

# values of the constants in the eq's
Lamb = 0
A = 1e-7
r_b = 5e10
n = 2
m = 2
H = 1
G = 6.67e-11  # Nm**2/kg**2 = m**3/kg*s**2
rho_i = 1  # 1e-26 #kg/m^3
args_list = (G, rho_i, r_b, n, m, A, H, Lamb)

# Initial value conditions
r_init = 1
R_init = 1
E_init = 1
M_init = 1
dMdr_init = 1
dRdr_init = 1
dEdr_init = 1
init_condt = [r_init, R_init, E_init, M_init, dMdr_init, dRdr_init, dEdr_init]