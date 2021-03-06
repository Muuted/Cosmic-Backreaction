# -------- Einstein_de_sitter_functions.py ---------------------
import numpy as np
from The_constants import *

def a_FLRW_lim(t,t_0):
    a_func = (t/t_0)**(2/3)
    return a_func

def Einstein_de_sitter(num_of_steps):#time_vec):
  
    def func_a_Ein_Sitter(a,G,rho_c0):

        const = np.sqrt(8*np.pi*G/(3))

        dadt = const*np.sqrt(rho_c0/a)

        return dadt


    def dSdt_a_ES(S,t,G,t_0,rho_c0):
        a = S

        dadt = func_a_Ein_Sitter(a, G,rho_c0)

        a_dot = (2/3)/(t_0**(2/3)*t**(1/3))
        a_func = (t/t_0)**(2/3)

        return dadt

    def func_rho_Ein_Sitter(a):
        Lamb, A, r_b, n, m, H_0, H_i, G, rho_c0, rho_i0, a_i, t_i, t_0,c= func_constants()

        rho_of_t = rho_c0/a**3

        return rho_of_t

    # ----------------------------- importing the constants--------------------------------
    Lamb, A, r_b, n, m, H_0, H_i, G, rho_c0, rho_i0, a_i, t_i, t_0,c= func_constants()

    init_cond = a_i
    num_of_steps = num_of_steps # given as argument for the function 
    time_vec = np.linspace(t_i,t_0, num_of_steps)

    ans_a_ES = scipy.integrate.odeint(func=dSdt_a_ES,
                                        y0=init_cond, t=time_vec, args=(G,t_0,rho_c0)
                                        )
    
    ans_a_ES = ans_a_ES.T

    a_ES = ans_a_ES[0]
    
    # finding rho
    rho_EdS = []#np.zeros(len(a_ES))
    for i in range(0,len(a_ES)):
        a_de_sitter = a_ES[i]
        #a_de_sitter = a_FLRW_lim(time_vec[i],t_0)
        rho_EdS.append(func_rho_Ein_Sitter(a_de_sitter))

    list = [a_ES, rho_EdS, time_vec]

    return list
"""
Lamb, A, r_b, n, m, H_0, H_i, G, rho_c0, rho_i0, a_i, t_i, t_0,c= func_constants()
num_of_steps = 10000

lit = Einstein_de_sitter(num_of_steps=num_of_steps)

a_ES = lit[0]
rho_ES = lit[1]
time_vec = lit[2]


print('time shape=',time_vec.shape)
print('a Es shape',a_ES.shape)
print('rho_max/rho_min=',max(rho_ES)/min(rho_ES))
print('a_max/a_min=',max(a_ES)/min(a_ES))


plt.figure()
plt.plot(time_vec,a_FLRW_lim(time_vec, t_0),label='analytic a')
plt.plot(time_vec,a_ES,'--',label='numeric a')
plt.legend()
plt.title('Einstein de Sitter model function')
#plt.show()


"""

def func_rho_Ein_Sitter(a):
    Lamb, A, r_b, n, m, H_0, H_i, G, rho_c0, rho_i0, a_i, t_i, t_0,c= func_constants()
    const = rho_c0

    rho_of_t = const/a**3

    return rho_of_t
