from LTB_model_functions import *
from Einstein_de_sitter_functions import *
#import scipy
#from scipy.integrate import solve_ivp
from The_constants import *
import plotly as py
import plotly.graph_objects as go

Lamb, A, r_b, n, m, H_0, H_i, G, rho_c0, rho_i0, a_i, t_i, t_0,c= func_constants()

num_steps = 10000 # number of steps between t_start and t_end
num_interations = 1 #number of r's
dt = 1e3        # time step
r = r_b*5         # distance
dr = 300        # change of r

EE  = 0
dEEdr  = 0
# Initial condition for the functions
RR = r*a_i
dRRdr = a_i
MM = (4*np.pi*G)/(3*c**2)*rho_i0*r**3*a_i**3
# --------------------- hvorfor skal det var 1/c^2 og ikke 1/c^4 
dMMdr = func_dMMdr(r, RR, 0, MM, 0, dRRdr, 0, G, rho_c0, r_b, n, m, A, H_i, Lamb, c)#- (4*np.pi*G)/(3*c**2)*rho_c0*3*r**2 # because its  dM/dr for the initial conditions 

EE = func_E(r, RR, 0, MM, dMMdr, dRRdr, 0, G, rho_c0, r_b, n, m, A, H_0, Lamb, c)
dEEdr = func_dEdr(r, RR, EE, MM, dMMdr, dRRdr, 0, G, rho_c0, r_b, n, m, A, H_0, Lamb, c)


# A list with all the arguments that is need to feed the functions.
args_list =[r, EE, dEEdr,MM, dMMdr, G, rho_c0, r_b, n, m, A, H_0, Lamb,c]

# Our time vector for the integration
time_tot = np.linspace(t_i,t_0,num_steps)


for i in range(0,num_interations):
    ans = []
    #The initial conditions are found for each r, and used in the ODE int integration
    init_cond_dRdt = [RR, dRRdr, rho_c0]#, MM ]

    ans_odeint = scipy.integrate.odeint(func_LTB_dSdt, 
            y0=init_cond_dRdt, 
            t=time_tot,
            args=(args_list,)
            )

    """
    ans_ivp = solve_ivp(
            dSdt_dRdrdt_ivp,t_span=(t_i,t_0),
            y0=init_cond_dRdt ,args=args_for_ivp, method ='RK45'
        )
    """



#FLRW_ans = func_FLRW_R_dRdr(time_tot, r, H_0)



# The results from our odeint  above
ans_odeint = ans_odeint.T

ans_RR = ans_odeint[0]
ans_dRdr = ans_odeint[1]
ans_rho = ans_odeint[2]


# Results for the Einstein de Sitter model 
a_ES, rho, rho_ES, time_vec = Einstein_de_sitter(num_of_steps=num_steps)
ans_a_ES = rho_ES

plt.figure()
plt.subplot(2,3,1)
plt.plot(time_tot,ans_RR/r, label='R(t,r)/r')
plt.xlabel('Gyr')
#plt.title('R/r')
plt.legend()

plt.subplot(2,3,2)
plt.plot(time_vec,a_ES,label='a(t) EdS')
plt.xlabel('Gyr')
plt.legend()
#plt.title('a(t) Ein de Sitter')

plt.subplot(2,3,3)
plt.plot(time_tot,ans_RR/r,label='R(t,r)/r')
plt.plot(time_vec,a_ES,'--',label='a(t) EdS')
plt.xlabel('Gyr')
#plt.title('R vs a EdS')
plt.legend()

plt.subplot(2,3,4)
plt.plot(time_tot,ans_dRdr,label='dRdr')
plt.xlabel('Gyr')
plt.legend()
#plt.title('dRdr')


plt.subplot(2,3,5)
#plt.plot(time_tot,ans_M,label='M(r)')
#plt.xlabel('Gyr')
#plt.legend()
#plt.title('M')


plt.subplot(2,3,6)
plt.plot(time_tot,ans_rho,label=r'$\rho$(t,r)')
plt.xlabel('Gyr')
plt.legend()
#plt.title('rho')

plt.show()
