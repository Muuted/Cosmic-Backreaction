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


# Initial condition for the functions
RR = r*a_i
dRRdr = a_i
MM = (4*np.pi*G)/(3*c**2)*rho_i0*r**3*a_i**3
# --------------------- hvorfor skal det var 1/c^2 og ikke 1/c^4 
dMMdr = - (4*np.pi*G)/(3*c**2)*rho_c0*3*r**2 # because its  dM/dr for the initial conditions 

EE = func_E(r, RR, 0, MM, dMMdr, dRRdr, 0, G, rho_c0, r_b, n, m, A, H_0, Lamb, c)
dEEdr = func_dEdr(r, RR, EE, MM, dMMdr, dRRdr, 0, G, rho_c0, r_b, n, m, A, H_0, Lamb, c)


# A list with all the arguments that is need to feed the functions.
args_list =[r, EE, dEEdr, dMMdr, G, rho_c0, r_b, n, m, A, H_0, Lamb,c]

# Our time vector for the integration
time_tot = np.linspace(t_i,t_0,num_steps)


for i in range(0,num_interations):
    ans = []
    #The initial conditions are found for each r, and used in the ODE int integration
    init_cond_dRdt = [RR, dRRdr, MM, rho_c0 ]

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
ans_M = ans_odeint[2]
ans_rho = ans_odeint[3]


# Results for the Einstein de Sitter model 
a_ES, rho, rho_ES, time_vec = Einstein_de_sitter(num_of_steps=num_steps)
ans_a_ES = rho_ES

plt.figure()
plt.subplot(2,3,1)
plt.plot(time_tot,ans_RR/r, label='R(t,r)/r')
#plt.title('R/r')
plt.legend()

plt.subplot(2,3,2)
plt.plot(time_vec,a_ES,label='a(t) EdS')
plt.legend()
#plt.title('a(t) Ein de Sitter')

plt.subplot(2,3,3)
plt.plot(time_tot,ans_RR/r,label='R(t,r)/r')
plt.plot(time_vec,a_ES,'--',label='a(t) EdS')
#plt.title('R vs a EdS')
plt.legend()

plt.subplot(2,3,4)
plt.plot(time_tot,ans_dRdr,label='dRdr')
plt.legend()
#plt.title('dRdr')


plt.subplot(2,3,5)
plt.plot(time_tot,ans_M,label='Mass')
plt.legend()
#plt.title('M')


plt.subplot(2,3,6)
plt.plot(time_tot,ans_rho,label='mass density')
plt.legend()
#plt.title('rho')


plt.show()
"""

fig_a = go.Figure()
fig_a.add_trace(
    go.Scatter(
            x=time_vec, y=a_ES,
            name='a(t)Einstein'
    )
)
fig_a.add_trace(
    go.Scatter(
        x=time_tot,
        y=ans_RR,
        name = "LTB R/r"
    )
)
fig_a.update_layout( title="a(t)",
        xaxis={'title': 'time [Gy]'},
        yaxis={
            'title': 'a(t)', 
            'tickformat':'e'
        },
        legend_title="Legend Title"
        )
fig_a.show()

"""
'''
fig_rho = go.Figure()
fig_rho.add_trace(
    go.Scatter(
    x=time_vec, y=rho
    )
)

fig_rho.update_layout( 
    title="rho",
    xaxis={'title': 'time [Gy]'},
    yaxis={'title': 'rho(t)',
            'tickformat':'e'
            }
)

fig_rho.show()
'''













"""
rho_list = []
rho_list_v2 = []
for i in range(0,len(R_vec)):

    RR = R_vec[i]
    dRRdr = dRdr_vec[i]
    args_list_ODE =[r, RR, EE, MM, dMMdr, dRRdr, dEEdr, G, rho_c0, r_b, n, m, A, H_0, Lamb,c]

    rho_list.append(
        func_rho( args_list_ODE )/rho_c0
    )

    EE = 0
    dEEdr = 0
    args_list_ODE_v2 =[r, RR, EE, MM, dMMdr, dRRdr, dEEdr, G, rho_c0, r_b, n, m, A, H_0, Lamb,c]

    rho_list_v2.append(
        func_rho( args_list_ODE_v2 )/rho_c0
    )
"""

def sdfsdf():
    plt.figure()
    plt.plot(time_tot,rho_list,label=r'$\rho (t) / \rho_{FLRW}$')
    plt.plot(time_tot,rho_list_v2,'--',label=r'$\rho (t) / \rho_{FLRW}$ _ v2')
    plt.xlim(-0.00001, 0.0050)
    plt.ylim(-0.1e-83 ,1.5e-83)
    plt.xlabel('Gy')
    plt.ylabel(r'$\rho$')
    plt.legend()
    """
    plt.figure()
    plt.plot(time_tot,rho_list_v2,label=r'$\rho (t) / \rho_{FLRW}$ _ v2')
    plt.xlim(-0.00001, 0.0050)
    plt.ylim(-0.1e-83 ,1.5e-83)
    plt.legend()
    """

    """
    plt.figure()
    plt.subplot(2,2,1)
    plt.plot(time_tot.T,R_vec,'k',label=r'$R(t,r)_{odeint}$')
    plt.plot(ans_ivp.t,ans_ivp.y[0],'--r',label=r'$R(t,r)_{ivp}$')
    plt.xlabel('Gy')
    plt.ylabel('R')
    plt.legend()

    plt.subplot(2,2,2)
    plt.plot(time_tot.T,dRdr_vec,'k',label=r'$\left(\frac{\partial R}{\partial r}\right)_{odeint}$')
    plt.plot(ans_ivp.t,ans_ivp.y[1],'--r',label=r'$ \left(\frac{\partial R}{\partial r} \right)_{ivp}$')
    plt.xlabel('Gy')
    plt.ylabel('dRdr')
    plt.legend()

    plt.subplot(2,2,3)
    plt.plot(time_tot, FLRW_ans[0],label=r'$R_{FLRW}$')
    plt.xlabel('Gy')
    plt.ylabel('R')
    plt.legend()

    plt.subplot(2,2,4)
    plt.plot(time_tot, FLRW_ans[1],label=r'$a_{FLRW}$')
    plt.xlabel('Gy')
    plt.ylabel('a - dRdr')
    plt.legend()

    plt.show()"""
