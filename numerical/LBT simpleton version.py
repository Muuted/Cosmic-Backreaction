from LTB_model_functions import *
from Einstein_de_sitter_functions import *
#import scipy
#from scipy.integrate import solve_ivp



# units for conversion
one_Mpc = 3.086e22 # m
one_Gy = 3.156e16 # s
one_year_in_sec = 60*60*24*365 # s
one_Gy = 3.156e16/one_year_in_sec # years
one_solar_mass = 1.989e30 #kg


# values of the constants in the eq's
Lamb = 0
A = 1e-7
r_b = 5e10
n = 2
m = 2

H_0 = 68 # km/s/Mpc -> Mpc/Gyr
H_0 = 68*1e3*one_Gy/one_Mpc # Mpc/Gy

G = 6.67e-11  #  m^3/kg/s^2
G = 6.67e-11*one_solar_mass*one_Gy**2/(one_Mpc**3) # Mpc^3/M_o*Gy^2

rho_FLRW =8.7e27 # kg/m^3 # This is the critical density Ryden in the back
rho_FLRW = 8.7e27*one_Mpc**3/one_solar_mass # M_o/Mpc^3


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
args_list =[r, RR, EE, 1, dMMdr, dRRdr, dEEdr, G, rho_FLRW, r_b, n, m, A, H_0, Lamb]
MM = 1#(4*np.pi/3)*r**3*rho_FLRW*(1-2*EE/(5*H_0**2))

# A list with all the arguments that is need to feed the functions.

args_list =[r, RR, EE, 1, dMMdr, dRRdr, dEEdr, G, rho_FLRW, r_b, n, m, A, H_0, Lamb]
time_tot = np.linspace(t_start,t_end,num_steps)#round(t_end/dt*10))


for i in range(0,num_interations):
    ans = []
    #The initial conditions are found for each r, and used in the ODE int integration
    EE = E(args_list)
    MM = M(args_list)
    dEEdr = dEdr(args_list)
    dMMdr= dMdr(args_list)

    # The constants under integration
    args_list_ODE =[r, RR, EE, MM, dMMdr, dRRdr, dEEdr, G, rho_FLRW, r_b, n, m, A, H_0, Lamb]

    args_for_ivp = (r, RR, EE, MM, dMMdr, dRRdr, dEEdr, G, rho_FLRW, r_b, n, m, A, H_0, Lamb)
    
    # the initial value for the function(s) that are being integrated
    init_cond_dRdt = [a_i*r, a_i]

    ans_odeint = scipy.integrate.odeint(
            dSdt_dRdrdt,t=time_tot,
            y0=init_cond_dRdt,args=(args_list_ODE,)
       )

    ans_ivp = solve_ivp(
            dSdt_dRdrdt_ivp,t_span=(t_start,t_end),
            y0=init_cond_dRdt ,args=args_for_ivp,method ='RK45'
        )



FLRW_ans = func_FLRW_R_dRdr(time_tot, r, H_0)


ans_odeint = ans_odeint.T
R_vec = ans_odeint[0]
dRdr_vec = ans_odeint[1]

print('H=',H_0, '\n',
    'G=',G,'\n',
    'rho_FLRW =',rho_FLRW,'\n',
    'min/max of time',min(time_tot),max(time_tot)
    )

rho_list = []
rho_list_v2 = []
for i in range(0,len(R_vec)):

    RR = R_vec[i]
    dRRdr = dRdr_vec[i]
    args_list_ODE =[r, RR, EE, MM, dMMdr, dRRdr, dEEdr, G, rho_FLRW, r_b, n, m, A, H_0, Lamb]

    rho_list.append(
        func_rho( args_list_ODE )/rho_FLRW
    )

    EE = 0
    dEEdr = 0
    args_list_ODE_v2 =[r, RR, EE, MM, dMMdr, dRRdr, dEEdr, G, rho_FLRW, r_b, n, m, A, H_0, Lamb]

    rho_list_v2.append(
        func_rho( args_list_ODE_v2 )/rho_FLRW
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

time_vec, a_ES, rho, rho_ES = Einstein_de_sitter()
ans_a_ES = rho_ES
fig_a = go.Figure()

fig_a.add_trace(
    go.Scatter(
            x=time_vec, y=ans_a_ES,#/max(ans_a_ES),
            name="a(t) Einstein"
    )
)
'''fig_a.add_trace(
    go.Scatter(
        x=time_tot,
        y=R_vec/max(R_vec),
        name = "LTB R/r"
    )
)'''
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