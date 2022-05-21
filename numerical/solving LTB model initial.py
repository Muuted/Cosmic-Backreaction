from logging import PlaceHolder
from LTB_model_functions import *
from Einstein_de_sitter_functions import *
from The_constants import *
import plotly as py
import plotly.graph_objects as go

Lamb, A, r_b, n, m, H_0, H_i, G, rho_c0, rho_i0, a_i, t_i, t_0,c= func_constants()

dt = 1e3        # time step
dr = 1        # change of r
r_i = 0.1         # distance, we start at origo 
num_steps = 10000 # number of steps between t_start and t_end
num_iterations = int(r_b/dr) #number of r's

# Our time vector for the integration
time_tot = np.linspace(t_i,t_0,num_steps)
radi_vec = np.linspace(r_b,r_i,num_iterations)

# Creating the lists for the data at different r values
ans_RR = [[] for x in range(num_iterations)]
ans_dRdr = [[] for x in range(num_iterations)]
ans_rho = [[] for x in range(num_iterations)]

Breakies = False

# increasing r for each iteration
r = 1 #radi_vec[i]

# Initial condition for the functions
RR = r*a_i
dRRdr = a_i

    
EE = func_E(r, r_b, n, m, A)
dEEdr = func_dEdr(r, r_b, n, m, A)
MM = func_M(r, EE, G, rho_c0, a_i, H_i, c)
dMMdr = func_dMMdr(r, EE, dEEdr, rho_c0, H_i, a_i, c) 

dRdt = func_dRdt(r, RR, EE, MM, G, c)
dRdrdt = func_dRdrdt(r, RR, EE, MM, dMMdr, dRRdr, dEEdr, G, c)
print('inner loop',dRdt,dRdrdt)
