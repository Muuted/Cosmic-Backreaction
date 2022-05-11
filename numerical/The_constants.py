import scipy
import scipy.integrate
import matplotlib.pyplot as plt
from scipy.integrate import quad
from scipy.integrate import solve_ivp
import numpy as np

def func_constants():
    
    # values of the constants in the eq's
    Lamb = 0
    A = 1e7
    r_b = 40
    n = 2
    m = 2

    # units for conversion
    m_pr_Mpc = 3.086e22 # m/Mpc
    one_year_in_sec = 60*60*24*365 # s
    s_pr_Gy = 3.156e16 # s/Gyr
    y_pr_Gy = 3.156e16/one_year_in_sec # years/Gyr
    kg_pr_solar_mass = 1.989e30 #kg/M_o
   
    
    # constants in normal units 
    H_0 = 68 # km/s/Mpc -> Mpc/Gyr
    G = 6.673e-11  #  m^3/kg/s^2
    c = 3e8 # m/s speed of light
    

    # constants in M_o and Mpc and Gyr
    H_0 = H_0*1e3*s_pr_Gy/m_pr_Mpc # 1/Gyr
    #H_0 = 0.07158985 # Hubble sonstant 70 km/s/Mpc in unit 1/Gyr
    G = G*kg_pr_solar_mass*s_pr_Gy**2/(m_pr_Mpc**3) # Mpc^3/M_o*Gy^2
    c = c*s_pr_Gy/m_pr_Mpc # Mpc/Gyr


    a_i = 1/1100    #initial scale factor.
    t_0 = 2/(3*H_0) # our time in Gyr
    t_i = t_0 * a_i**(3/2) #start time in Gyr

    rho_c0 = 3*H_0**2/(8*np.pi*G*a_i**3) # M_0 / Mpc^3 should be the right units as G and H_0 are  
    #rho_c0 = 1.28e11 # Ryden page : 57 eq (4.32)

    
    H_i = np.sqrt(8*np.pi*G*rho_c0/3) # initial

    
    const_list = [Lamb, A, r_b, n, m, H_0, G, rho_c0, a_i, t_i, t_0,c]
    return const_list


Lamb, A, r_b, n, m, H_0, G, rho_c0, a_i, t_i, t_0, c= func_constants()
