import scipy
import scipy.integrate

import matplotlib.pyplot as plt

from scipy.integrate import quad

from scipy.integrate import solve_ivp

import numpy as np

#from constants import *

# units for conversion
one_Mpc = 3.086e22 # m
one_Gy = 3.156e16 # s
one_Gy = 3.156e16/60*60*24*365 # years
one_solar_mass = 1.989e30 #kg



def func_a_Ein_Sitter(a,G):

    const = 8*np.pi*G/3
    dadt = np.sqrt(
        const/a
    )

    return dadt

def func_rho_Ein_Sitter(a):

    const = 1

    drhodt = const/a**3

    return drhodt