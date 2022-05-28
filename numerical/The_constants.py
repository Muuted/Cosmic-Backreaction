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
    y_pr_Gy = 3.156e16/one_year_in_sec # years/Gyr
    kg_pr_solar_mass = 1.989e30 #kg/M_o


    # values of the constants in the eq's
    Lamb = 0
    A = 1.3e-7  # [A] = 1/Mpc^2
    # Wrong A below
    #A = 1.3e-14
    
    r_b = 50
    n = 2
    m = 2

    
    # constants in normal units 
    H_0 = 70 # km/s/Mpc -> Mpc/Gyr
    G = 6.673e-11  #  m^3/kg/s^2
    c = 3e8#299792458 #3e8 # m/s speed of light
    #H_0 = H_0*1e3*(1/m_pr_Mpc)    #Si 1/s

    # constants in M_o and Mpc and Gyr
    H_0 = H_0*1e3*s_pr_Gy/m_pr_Mpc # 1/Gyr
    G = G*kg_pr_solar_mass*s_pr_Gy**2/(m_pr_Mpc**3) # Mpc^3/M_o*Gy^2
    c = c*s_pr_Gy/m_pr_Mpc # Mpc/Gyr


    a_i = 1/1100    #initial scale factor.
    t_0 = 2/(3*H_0) # our time in Gyr
    t_i = t_0 * a_i**(3/2) #start time in Gyr

    rho_c0 = 3*H_0**2/(8*np.pi*G) # M_0 / Mpc^3 should be the right units as G and H_0 are  
    #rho_c0 = 3*H_0**2/(8*np.pi*G)*c**2 # M_0 / Mpc^3 should be the right units as G and H_0 are  
    
    rho_i0 = rho_c0/a_i**3
    H_i = np.sqrt(8*np.pi*G*rho_i0/3) # initial

    # Wrong H_i
    #H_i = np.sqrt(8*np.pi*G*rho_c0/3) # initial

      
    const_list = [Lamb, A, r_b, n, m, H_0, H_i, G, rho_c0, rho_i0, a_i, t_i, t_0,c]
    return const_list


Lamb, A, r_b, n, m, H_0, H_i, G, rho_c0, rho_i0, a_i, t_i, t_0,c= func_constants()
#print('t_i=',t_i/(60*60*24*365),'\n','t_0=',t_0/(60*60*24*365))
#print(c,'\n',G*1e15,'\n',H_0,'\n',a_i,'\n',t_i,'\n',rho_c0/1e12)

print(t_0)
#print(H_i,H_0)

H0 = 0.07158985 # Hubble sonstant 70 km/s/Mpc in unit 1/Gyr
t0 = 2./3/H0 # assume EdS background
c = 306.60139 # speed of light in vacuum in units Mpc/Gyr
M_sun = 1.98855e30
Gyr = 3.15576e16
Mpc = 3.0856778e22
r_b = 40
k_max = 1.3e-7
a_i = 1.0/1100
G = 6.673e-11*M_sun*Gyr*Gyr/Mpc/Mpc/Mpc
rho_ic = 3.0*H0*H0/(8.0*np.pi*G)/a_i**3
Hi = ( 8*np.pi*G*rho_ic/3.0 )**(1/2)

#print(Hi,H0)