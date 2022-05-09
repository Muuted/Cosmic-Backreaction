import scipy
import scipy.integrate

import matplotlib.pyplot as plt

from scipy.integrate import quad

from scipy.integrate import solve_ivp

import numpy as np

from constants import *

import plotly as py
import plotly.graph_objects as go

def Einstein_de_sitter(time_vec):
  
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
    a_i = 1# 1/1100
    rho_i = 3*H_0**2/(8*np.pi*G)

    init_cond = [a_i, rho_i]

    t_0 = 14e9 # years
    t_0 = t_0/one_Gy # Gigayears
    t_i = t_0 * a_i**(3/2)

    num_of_steps = 100
    #time_vec = np.linspace(t_i,t_0, num_of_steps)

    ans_a_ES = scipy.integrate.odeint(func=dSdt_a_ES,
                                        y0=init_cond, t=time_vec, args=(G,)
                                        )
    
    ans_a_ES = ans_a_ES.T

    a_ES = ans_a_ES[0]
    rho_ES = ans_a_ES[1]

    # finding rho
    rho = np.zeros(len(a_ES))
    #rho = []
    for i in range(len(a_ES)):
        a_de_sitter = a_ES[i]
        rho[i]=func_rho_Ein_Sitter(a_de_sitter)

    print('time shape=',time_vec.shape)
    list = [a_ES, rho, rho_ES]
    return list
"""
lit = Einstein_de_sitter()


a_ES = lit[0]
rho = lit[1]
rho_ES = lit[2]
t_i = lit[3]
t_0 = lit[4]
num_of_steps = lit[5]
t_0 = 14e9 # years
t_0 = t_0/one_Gy # Gigayears
t_i = t_0 * a_i**(3/2)

num_of_steps = 100
time_vec = np.linspace(t_i,t_0, num_of_steps)

print('time shape=',time_vec.shape)
print('a Es shape',a_ES.shape)
print('rho_max/rho_min=',max(rho_ES)/min(rho_ES))
print('a_max/a_min=',max(a_ES)/min(a_ES))

plt.figure()
plt.plot(time_vec,rho_ES)
plt.plot(time_vec,a_ES)
plt.show()"""