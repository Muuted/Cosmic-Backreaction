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
    MM = 4 * np.pi * RR ** 3 * rho_i / 3 * (1 - 2 * EE / (5 * H ** 2 ** RR ** 2))
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


def dSdt(r, RR, EE, MM, dMMdr, dRRdr, dEEdr, G, rho_i, r_b, n, m, A, H, Lamb):
    sol_EE = E(r, RR, EE, MM, dMMdr, dRRdr, dEEdr, G, rho_i, r_b, n, m, A, H, Lamb)
    sol_MM = M(r, RR, EE, MM, dMMdr, dRRdr, dEEdr, G, rho_i, r_b, n, m, A, H, Lamb)
    sol_dEEdr = dEdr(r, RR, EE, MM, dMMdr, dRRdr, dEEdr, G, rho_i, r_b, n, m, A, H, Lamb)
    sol_dMMdr = dMdr(r, RR, EE, MM, dMMdr, dRRdr, dEEdr, G, rho_i, r_b, n, m, A, H, Lamb)
    sol_dRdt = dRdt(r, RR, EE, MM, dMMdr, dRRdr, dEEdr, G, rho_i, r_b, n, m, A, H, Lamb)
    sol_rho = rho(r, RR, EE, MM, dMMdr, dRRdr, dEEdr, G, rho_i, r_b, n, m, A, H, Lamb)
    return [sol_EE, sol_MM, sol_dEEdr, sol_dMMdr, sol_dRdt, sol_rho]


def dSdt2(r, RR, EE, MM, dMMdr, dRRdr, dEEdr, G, rho_i, r_b, n, m, A, H, Lamb):
    sol_EE = E(r, RR, EE, MM, dMMdr, dRRdr, dEEdr, G, rho_i, r_b, n, m, A, H, Lamb)
    sol_MM = M(r, RR, EE, MM, dMMdr, dRRdr, dEEdr, G, rho_i, r_b, n, m, A, H, Lamb)
    sol_dRRdt = dRdt(r, RR, EE, MM, dMMdr, dRRdr, dEEdr, G, rho_i, r_b, n, m, A, H, Lamb)


Lamb = 0
A = 1e-7
r_b = 1e6
n = 2
m = 2
H = 1
G = 6.67 * 1e27
rho_i = 1
sol = solve_ivp(dSdt, t_span=[0, 100], y0=[1, 1, 1, 1, 1, 1, 1],
                args=(G, rho_i, r_b, n, m, A, H, Lamb,1, 1, 1, 1, 1), method='RK45')
