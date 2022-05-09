import scipy
import scipy.integrate

import matplotlib.pyplot as plt

from scipy.integrate import quad

from scipy.integrate import solve_ivp

import numpy as np

#from constants import *

from constants import *

import plotly as py
import plotly.graph_objects as go

def Einstein_de_sitter():
    # units for conversion
    one_Mpc = 3.086e22 # m
    one_Gy = 3.156e16 # s
    one_Gy = one_Gy/(60*60*24*365) # years
    print('Gy in years =',one_Gy)
    one_solar_mass = 1.989e30 #kg


    def func_a_Ein_Sitter(a,G):

        const = np.sqrt(8*np.pi*G/3)

        dadt = const*np.sqrt(1/a)

        return dadt
    
    def func_drhodt(a,dadt):

        drhodt = -dadt/(a**4)

        return drhodt

    def dSdt_a_ES(S,t,G):
        a,rho = S

        dadt = func_a_Ein_Sitter(a, G)

        drhodt = func_drhodt(a,dadt)

        return [dadt, drhodt]

    def func_rho_Ein_Sitter(a):

        const = 1

        rho_of_t = const/a**3

        return rho_of_t


    # initial conditions 
    a_i = 1/1100
    rho_i = 3*H_0**2/(8*np.pi*G)

    init_cond = [a_i, rho_i]
    t_0 = 14e9 # years
    t_0 = t_0/one_Gy # Gigayears

    t_i = t_0 * a_i**(3/2)

    
    num_of_steps = 10000
    time_vec = np.linspace(t_i,t_0, num_of_steps)

    ans_a_ES = scipy.integrate.odeint(func=dSdt_a_ES,
                                        y0=init_cond, t=time_vec, args=(G,)
                                        )


    
    ans_a_ES = ans_a_ES.T

    a_ES = ans_a_ES[0].T
    rho_ES = ans_a_ES[1].T
    


    # finding rho
    rho = []
    for i in range(len(ans_a_ES)):
        a_de_sitter = a_ES[i]
        rho.append(func_rho_Ein_Sitter(a_de_sitter))


    print(max(rho_ES)/min(rho_ES))
    print(max(a_ES)/min(a_ES))

    plt.figure()
    plt.plot(time_vec,rho_ES)
    plt.show()
    list = [time_vec, a_ES, rho, rho_ES]
    return list

time_vec, ans_a_ES, rho = Einstein_de_sitter()

