from LTB_model_functions import *
from Einstein_de_sitter_functions import *
from The_constants import *
import plotly as py
import plotly.graph_objects as go

Lamb, A, r_b, n, m, H_0, H_i, G, rho_c0, rho_i0, a_i, t_i, t_0,c= func_constants()

num_steps = 10000 # number of steps between t_start and t_end
num_interations = 1 #number of r's
dt = 1e3        # time step
r = r_b*5         # distance
dr = 300        # change of r

# Initial condition for the functions
RR = r*a_i
dRRdr = a_i

EE = func_E(r, r_b, n, m, A)
dEEdr = func_dEdr(r, r_b, n, m, A)
MM = func_M(r, EE, G, rho_c0, a_i, H_i, c)#(4*np.pi*G)/(3*c**2)*rho_i0*r**3*a_i**3
# --------------------- hvorfor skal det var 1/c^2 og ikke 1/c^4 
dMMdr = func_dMMdr(r, EE, dEEdr, rho_c0, H_i, a_i, c) #- (4*np.pi*G)/(3*c**2)*rho_c0*3*r**2 # because its  dM/dr for the initial conditions 


# A list with all the arguments that is need to feed the functions.
args_list =[r, EE, dEEdr,MM, dMMdr, G, c]
# Our time vector for the integration
time_tot = np.linspace(t_i,t_0,num_steps)
# The initial conditions are found for each r, 
# and used in the ODE int integration
init_cond_dRdt = [RR, dRRdr]
# ----------------------------Integration start ------------------------------------------------------------
ans_odeint = scipy.integrate.odeint(func_LTB_dSdt, 
        y0=init_cond_dRdt, 
        t=time_tot,
        args=(args_list,)
        )
# The results from our odeint  above
ans_odeint = ans_odeint.T
# --------------------------Integration done ------------------------------------------------------------

ans_RR = ans_odeint[0]
ans_dRdr = ans_odeint[1]



# Results for the Einstein de Sitter model 
a_ES, rho, rho_ES, time_vec = Einstein_de_sitter(num_of_steps=num_steps)
ans_a_ES = rho_ES

plt.figure()
plt.subplot(2,3,1)
plt.plot(time_tot,ans_RR/r, label=r'$\dfrac{R(t,r)}{r}$')
plt.xlabel('Gyr')
#plt.title('R/r')
plt.legend()


plt.subplot(2,3,2)
plt.plot(time_vec,a_ES,label=r'$a_{EdS}$')
plt.xlabel('Gyr')
plt.legend()
#plt.title('a(t) Ein de Sitter')

plt.subplot(2,3,3)
plt.plot(time_tot,ans_RR/r,label=r'$\dfrac{R(t,r)}{r}$')
plt.plot(time_vec,a_ES,'--',label=r'$a_{EdS}$')
plt.xlabel('Gyr')
#plt.title('R vs a EdS')
plt.legend()

plt.subplot(2,3,4)
plt.plot(time_tot,ans_dRdr,label=r'$\dfrac{\partial R(t,r)}{\partial r}$')
plt.xlabel('Gyr')
plt.legend()
#plt.title('dRdr')


plt.subplot(2,3,5)
plt.plot(time_tot,func_rho(r, ans_RR, dMMdr, ans_dRdr, rho_c0),'-o',label=r'$\rho_{LTB}$(t,r)')
plt.plot(time_vec,rho,label=r'$\rho_{EdS}$')

plt.xlabel('Gyr')
plt.legend()
#plt.title('rho')


plt.show()
