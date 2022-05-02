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

def dSdt(S, t, G, rho_i, r_b, n, m, A, H, Lamb):
    r, RR, EE, MM, dMMdr, dRRdr, dEEdr = S

    sol_EE = E(r, RR, EE, MM, dMMdr, dRRdr, dEEdr, G, rho_i, r_b, n, m, A, H, Lamb)
    sol_MM = M(r, RR, EE, MM, dMMdr, dRRdr, dEEdr, G, rho_i, r_b, n, m, A, H, Lamb)
    sol_dEEdr = dEdr(r, RR, EE, MM, dMMdr, dRRdr, dEEdr, G, rho_i, r_b, n, m, A, H, Lamb)
    sol_dMMdr = dMdr(r, RR, EE, MM, dMMdr, dRRdr, dEEdr, G, rho_i, r_b, n, m, A, H, Lamb)
    sol_dRdt = dRdt(r, RR, EE, MM, dMMdr, dRRdr, dEEdr, G, rho_i, r_b, n, m, A, H, Lamb)
    sol_rho = rho(r, RR, EE, MM, dMMdr, dRRdr, dEEdr, G, rho_i, r_b, n, m, A, H, Lamb)
    return [r, sol_EE, sol_MM, sol_dEEdr, sol_dMMdr, sol_dRdt, sol_rho]

'''
def dSdt2(S, r, G, rho_i, r_b, n, m, A, H, Lamb):
    RR, EE, MM, = S

    sol_dEEdr = dEdr(r, RR, EE, MM, dMMdr, dRRdr, dEEdr, G, rho_i, r_b, n, m, A, H, Lamb)
    sol_dMMdr = dMdr(r, RR, EE, MM, dMMdr, dRRdr, dEEdr, G, rho_i, r_b, n, m, A, H, Lamb)
    
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

# sol = solve_ivp(dSdt, t_span=[0, 100], y0=init_condt, args=(args_list))

# using the odeint to find the solution to the equations.
time_tot = np.linspace(0, 10, 10)
ans = scipy.integrate.odeint(dSdt, y0=init_condt, t=time_tot, args=args_list)

print(ans.shape, ans.shape[0], ans.shape[1])
# transpose the array because we have that the a arrays
# are of the format [[ r0, sol0 ,...],[r1,sol1..],[r2, sol2, ...] ,... ]
# so we transpose to have ans[0] to be the full array of r values.
ans = ans.T
print(ans.shape, ans.shape[0], ans.shape[1])
print(ans)

plt.figure()
plt.plot(ans[0], label='r')
plt.plot(ans[1], label='E(r)')
plt.plot(ans[2], label='M(r)')
plt.plot(ans[3], label='dEdr')
#plt.plot(ans[4], label='dMdr')
plt.plot(ans[5], label='dRdt')
plt.plot(ans[6], label='rho')
plt.legend()

plt.show()
