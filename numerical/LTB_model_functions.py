# ----------------------------------- LTB_model_functions.py --------------------------
import numpy as np

def func_E(r, r_b, n, m, A):

    if r <= r_b:

        EE = A * r ** 2 * ((r / r_b) ** n - 1) ** m

    else:

        EE = 0

    return EE

def func_dEdr(r,  r_b, n, m, A):

    if r <= r_b :
        dEEdr1 = 2 * A * r * ((r / r_b) ** n - 1) ** m

        dEEdr2 = ((m* n * A )/(r_b**n)) * r ** (n + 1) * ((r / r_b) ** n - 1) ** (m - 1)

        dEEdr = dEEdr1 + dEEdr2
    else:
        dEEdr = 0
    return dEEdr

def func_rho(RR, dMMdr, dRRdr, c):

    kappa = (8*np.pi)
    rho = (2/kappa)*(dMMdr/(RR**2*dRRdr))

    return rho

def func_M(r, EE, G, rho_c0,a_i,H_i, c):
    
    C_1 = 4*np.pi*rho_c0/3
    C_2 = (6*c**2/(5*a_i**2*H_i**2))

    ana_M = C_1*r**3*(1 - C_2*EE/(r**2))

    return ana_M

def func_dMMdr(r, EE, dEEdr, rho_c0, H_i,a_i, c):
    
    C_1 = 4*np.pi*rho_c0
    C_2 = (6*c**2)/(5*a_i**2*H_i**2)

    dMMdr = C_1*(r**2 - C_2*(EE + r*dEEdr))

    return dMMdr

def func_dRdt(RR, EE, MM, G, c):
    
    dRRdt = c*np.sqrt(
            2 * G*MM /(c**2* RR) + 2 * EE
        )

    return dRRdt

def func_dRdrdt(RR, EE, MM, dMMdr, dRRdr, dEEdr, G, c):
    
    sqrts = np.sqrt(
           2*G*MM/(c**2*RR) +2*EE #+ (Lamb/3)*RR**2 
      )

    extra = (G/c)*(dMMdr/RR - MM*dRRdr/(RR**2)) + c*dEEdr #+ Lamb*RR*dRRdr/3

    return extra/sqrts

def func_LTB_dSdt(S,t,p):    
    RR, dRRdr = S
    
    EE, dEEdr, MM, dMMdr, G, c = p

    dRdt = func_dRdt(RR, EE, MM, G, c)
    dRdrdt =  func_dRdrdt(RR, EE, MM, dMMdr, dRRdr, dEEdr, G, c)
    
    return_list = [dRdt, dRdrdt]
    return  return_list
    