import scipy
import scipy.integrate

import matplotlib.pyplot as plt

from scipy.integrate import quad

from scipy.integrate import solve_ivp

import numpy as np


# units for conversion
one_Mpc = 3.086e22 # m
one_Gy = 3.156e16 # s
one_Gy = 3.156e16/60*60*24*365 # years
one_solar_mass = 1.989e30 #kg


# values of the constants in the eq's
Lamb = 0
A = 1e7
r_b = 5e10
n = 2
m = 2

H_0 = 68 # km/s/Mpc -> Mpc/Gyr
H_0 = 68*1e3*one_Gy#/one_Mpc # Mpc/Gy
print(H_0)
H_0 = 1/13.4 # 1/Gy
print(H_0)
G = 6.67e-11  #  m^3/kg/s^2
G = G*one_solar_mass*one_Gy**2/(one_Mpc**3) # Mpc^3/M_o*Gy^2

rho_FLRW = 3*H_0**2/(8*np.pi*G) #8.7e27 # kg/m^3 # This is the critical density Ryden in the back
print(G,rho_FLRW)
rho_FLRW = rho_FLRW*one_Mpc**3/one_solar_mass # M_o/Mpc^3
print(rho_FLRW)

# The unit for rho_FLRW is fond on the website: 
# http://astroweb.case.edu/ssm/astroconstants.html
#   #rho_FLRW = 1.5e-7
# right above the Galactic untis


a_i = 1/1100    #initial scale factor.
t_end = 14e9    # end time (our time) in years
t_end = t_end/one_Gy # end time in Gy

t_start= 370e3   # start time years  found on -> https://en.wikipedia.org/wiki/Chronology_of_the_universe
                 # under: The Early Universe
t_start = t_start/one_Gy # start time in Gy


# Normalizing the time vector.
# So is ends at one, our time
t_start = t_start/t_end
t_end = t_end/t_end

