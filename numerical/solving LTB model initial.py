from scipy.integrate import quad
import numpy as np


def E(r, r_b, A, n, m):
    if r <= r_b:
        EE = A * r ** 2 * ((r / r_b) ** n - 1) ** m
    else:
        EE = 0
    return EE


def M(E, R, rho_i, H):
    MM = 4 * np.pi * R ** 3 * rho_i / 3 * (1 - 2 * E / (5 * H ** 2 ** R ** 2))
    return MM


def dEdr(r, A, r_b, n, m):
    dEEdr1 = 2 * A * r * ((r / r_b) ** n - 1) ** m
    dEEdr2 = m * n * A * r ** (n + 1) * ((r / r_b) ** n - 1) ** (m - 1)
    return dEEdr1 + dEEdr2


def dMdr(RR, dRRdr, rho_i, E, H, dEEdr):
    dMMdr1 = 12 * np.pi * RR ** 2 * dRRdr * rho_i / 3
    dMMdr2 = -(8 * np.pi * rho_i / 15 * H ** 2) * (E * dRRdr + RR * dEEdr)
    return dMMdr1 + dMMdr2


def dRdt(RR, MM, EE, Lamb):
    dRRdt = np.sqrt(2 * MM / RR + 2 * EE + (Lamb / 3) * RR ** 2)
    return dRRdt


def rho(RR, dMMdr, dRRdr,G):
    rrho = (1 / (4 * np.pi * G)) * dMMdr / (RR ** 2 * dRRdr)
    return rrho


Lamb = 0
A = 1e-7
r_b = 1e6
n = 2
m = 2
H = 1


sols = quad()


