from LTB_model_functions import *
from Einstein_de_sitter_functions import *
#import scipy
#from scipy.integrate import solve_ivp
from The_constants import *
import plotly as py
import plotly.graph_objects as go

Lamb, A, r_b, n, m, H_0, G, rho_c0, a_i, t_i, t_0, c = func_constants()

num_steps = 10000 # number of steps between t_start and t_end
num_interations = 1 #number of r's
dt = 1e3        # time step
r = 0.1         # distance
dr = 300        # change of r


RR = r*a_i
dRRdr = a_i 
EE = 0
dMMdr = 1
dEEdr = 0
args_list =[r, RR, EE, 1, dMMdr, dRRdr, dEEdr, G, rho_c0, r_b, n, m, A, H_0, Lamb]
MM = 1#(4*np.pi/3)*r**3*rho_FLRW*(1-2*EE/(5*H_0**2))

# A list with all the arguments that is need to feed the functions.
args_list =[r, RR, EE, 1, dMMdr, dRRdr, dEEdr, G, rho_c0, r_b, n, m, A, H_0, Lamb]

# Our time vector for the integration
time_tot = np.linspace(t_i,t_0,num_steps)


for i in range(0,num_interations):
    ans = []
    #The initial conditions are found for each r, and used in the ODE int integration
    EE = E(args_list)
    MM = M(args_list)
    dEEdr = dEdr(args_list)
    dMMdr= dMdr(args_list)

    # The constants under integration
    args_list_ODE =[r, RR, EE, MM, dMMdr, dRRdr, dEEdr, G, rho_c0, r_b, n, m, A, H_0, Lamb]

    args_for_ivp = (r, RR, EE, MM, dMMdr, dRRdr, dEEdr, G, rho_c0, r_b, n, m, A, H_0, Lamb)
    
    # the initial value for the function(s) that are being integrated
    init_cond_dRdt = [a_i*r, a_i]

    ans_odeint = scipy.integrate.odeint(
            dSdt_dRdrdt,t=time_tot,
            y0=init_cond_dRdt,args=(args_list_ODE,)
       )

    ans_ivp = solve_ivp(
            dSdt_dRdrdt_ivp,t_span=(t_i,t_0),
            y0=init_cond_dRdt ,args=args_for_ivp, method ='RK45'
        )



FLRW_ans = func_FLRW_R_dRdr(time_tot, r, H_0)


ans_odeint = ans_odeint.T
R_vec = ans_odeint[0]
dRdr_vec = ans_odeint[1]

print('H=',H_0, '\n',
    'G=',G,'\n',
    'rho_FLRW =',rho_c0,'\n',
    'min/max of time',min(time_tot),max(time_tot)
    )

rho_list = []
rho_list_v2 = []
for i in range(0,len(R_vec)):

    RR = R_vec[i]
    dRRdr = dRdr_vec[i]
    args_list_ODE =[r, RR, EE, MM, dMMdr, dRRdr, dEEdr, G, rho_c0, r_b, n, m, A, H_0, Lamb]

    rho_list.append(
        func_rho( args_list_ODE )/rho_c0
    )

    EE = 0
    dEEdr = 0
    args_list_ODE_v2 =[r, RR, EE, MM, dMMdr, dRRdr, dEEdr, G, rho_c0, r_b, n, m, A, H_0, Lamb]

    rho_list_v2.append(
        func_rho( args_list_ODE_v2 )/rho_c0
    )



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

"""plt.figure()
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

a_ES, rho, rho_ES, time_vec = Einstein_de_sitter(num_of_steps=num_steps)
ans_a_ES = rho_ES
"""t_0 = 14e9 # years
t_0 = t_0/one_Gy # Gigayears
t_i = t_0 * a_i**(3/2)
num_of_steps = 100"""
#time_vec = np.linspace(t_i,t_0, num_of_steps)


fig_a = go.Figure()
fig_a.add_trace(
    go.Scatter(
            x=time_vec, y=a_ES,#/max(ans_a_ES),
            name="a(t)Einstein"
    )
)
fig_a.add_trace(
    go.Scatter(
        x=time_tot,
        y=R_vec/max(R_vec),
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