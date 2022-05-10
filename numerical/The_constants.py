import scipy
import scipy.integrate
import matplotlib.pyplot as plt
from scipy.integrate import quad
from scipy.integrate import solve_ivp
import numpy as np

def func_constants():
    
    # units for conversion
    m_pr_Mpc = 3.086e22 # m/Mpc
    one_year_in_sec = 60*60*24*365 # s
    s_pr_Gy = 3.156e16 # s/Gyr
    s_pr_Gy = 3.156e16/one_year_in_sec # years/Gyr
    kg_pr_solar_mass = 1.989e30 #kg/M_o

    # values of the constants in the eq's
    Lamb = 0
    A = 1e7
    r_b = 5e10
    n = 2
    m = 2

    # constants in normal units 
    H_0 = 68 # km/s/Mpc -> Mpc/Gyr
    G = 6.67e-11  #  m^3/kg/s^2

    # constants in M_o and Mpc and Gyr
    H_0 = 68*1e3*s_pr_Gy/m_pr_Mpc # 1/Gyr
    G = G*kg_pr_solar_mass*s_pr_Gy**2/(m_pr_Mpc**3) # Mpc^3/M_o*Gy^2
    rho_c0 = 3*H_0**2/(8*np.pi*G) # M_0 / Mpc^3 should be the right units as G and H_0 are  
    # rho is ~ 1.28e11 M_0 / Mpc^3 ,just like in Ryden.


    a_i = 1/1100    #initial scale factor.
    t_0 = 2/(3*H_0) # our time in Gyr
    t_i = t_0 * a_i**(3/2) #start time in Gyr

    
    const_list = [Lamb, A, r_b, n, m, H_0, G, rho_c0, a_i, t_i, t_0]
    return const_list


#Lamb, A, r_b, n, m, H_0, G, rho_c0, a_i, t_i, t_0 = constants()
#print(Lamb, A, r_b, n, m, H_0, G, rho_c0, a_i, t_i, t_0)