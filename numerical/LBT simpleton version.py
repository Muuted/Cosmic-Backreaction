from LTB_model_functions import *
#import scipy
#from scipy.integrate import solve_ivp

def dSdt_dRdrdt(S,t,p):
    R, dRdr = S
    #r, RR, EE, MM, dMMdr, dRRdr, dEEdr, G, rho_FLRW, r_b, n, m, A, H, Lamb = p
    
    ans_dRdr = ddRdrdt(r, RR, EE, MM, dMMdr, dRRdr, dEEdr, G, rho_FLRW, r_b, n, m, A, H, Lamb)
    ans_R = dRdt(r, RR, EE, MM, dMMdr, dRRdr, dEEdr, G, rho_FLRW, r_b, n, m, A, H, Lamb)
    
    return [ans_R, ans_dRdr]

def dSdt_dRdrdt_ivp(t,S,r, RR, EE, MM, dMMdr, dRRdr, dEEdr, G, rho_FLRW, r_b, n, m, A, H, Lamb):
    R, dRdr = S
    #r, RR, EE, MM, dMMdr, dRRdr, dEEdr, G, rho_FLRW, r_b, n, m, A, H, Lamb = p
    
    ans_dRdr = ddRdrdt(r, RR, EE, MM, dMMdr, dRRdr, dEEdr, G, rho_FLRW, r_b, n, m, A, H, Lamb)
    ans_R = dRdt(r, RR, EE, MM, dMMdr, dRRdr, dEEdr, G, rho_FLRW, r_b, n, m, A, H, Lamb)
    
    return [ans_R, ans_dRdr]

# units for conversion
one_Mpc = 3.086e22 # m
one_Gy = 3.156e16 #s
one_solar_mass = 1.989e30 #kg

# values of the constants in the eq's
Lamb = 0
A = 1e7
r_b = 5e10
n = 2
m = 2

H = 68 # km/s/Mpc -> Mpc/Gyr
H = 68*1e3*one_Gy/one_Mpc # Mpc/Gy

G = 6.67e-11  #  m^3/kg/s^2
G = 6.67e-11*one_solar_mass*one_Gy**2/(one_Mpc**3) # Mpc^3/M_o*Gy^2

rho_FLRW =8.7e27 # kg/m^3 # This is the critical density Ryden in the back
rho_FLRW = 1.5e-7 #M_o/Mpc^3 # 8.7e27*one_Mpc**3/one_solar_mass # M_o/Mpc^3
# The unit for rho_FLRW is fond on the website:
# http://astroweb.case.edu/ssm/astroconstants.html
# right above the Galactic untis

a_i = 1/1100    #initial scale factor.
t_start= 0    # start time
t_end = 14e6    # end time
num_steps = 10000 # number of steps between t_start and t_end
num_interations = 1 #number of r's
dt = 1e3        # time step
r = r_b         # distance
dr = 300        # change of r


dRRdr = a_i 
MM = 1
EE = 0
dMMdr = 1
dRRdr = 1
dEEdr = 1
RR = r*a_i


# A list with all the arguments that is need to feed the functions.
args_list =[r, RR, EE, MM, dMMdr, dRRdr, dEEdr, G, rho_FLRW, r_b, n, m, A, H, Lamb]


time_tot = np.linspace(t_start,t_end,num_steps)#round(t_end/dt*10))


ivp = False 

for i in range(0,num_interations):
    ans = []
    #The initial conditions are found for each r, and used in the ODE int integration
    
    EE = E(args_list)
    MM = M(args_list)

    dEEdr = dEdr(args_list)
    dMMdr= dMdr(args_list)

    # The constants under integration
    args_list_ODE =[r, RR, EE, MM, dMMdr, dRRdr, dEEdr, G, rho_FLRW, r_b, n, m, A, H, Lamb]

    args_for_ivp = (r, RR, EE, MM, dMMdr, dRRdr, dEEdr, G, rho_FLRW, r_b, n, m, A, H, Lamb)
    

    # the initial value for the function(s) that are being integrated
    init_cond_dRdt = [a_i*r, a_i]
    
    r += dr
    
    if ivp == False:
        ans = scipy.integrate.odeint(dSdt_dRdrdt,t=time_tot,y0=init_cond_dRdt,
            args=(args_list_ODE,)
            )
    if ivp == True:
        ans = solve_ivp(dSdt_dRdrdt_ivp,t_span=(t_start,t_end),
                        y0=init_cond_dRdt ,args=args_for_ivp,
                        method ='LSODA'
                        #,t_eval = time_tot
                        )
        

if ivp == False:
    ans = ans.T
    R_vec = ans[0]
    print('len(R)=',len(R_vec))
    dRdr_vec = ans[1]
    print('len(dRdr)=',len(dRdr_vec))
    print(R_vec)
    print('R_max /R_min=',max(R_vec)/min(R_vec), '\n',
        'dRdr_max/dRdr_min =',max(dRdr_vec)/min(dRdr_vec))
    print('H=',H, '\n',
        'G=',G,'\n',
        'rho_FLRW =',rho_FLRW
        )
    plt.figure()
    plt.plot(time_tot.T,R_vec,'-',label='R(t,r)')
    plt.legend()

    plt.figure()
    plt.plot(time_tot.T,dRdr_vec,'-',label=r'$\frac{\partial R}{\partial r}$')
    plt.legend()
    

if ivp == True:

    plt.figure()

    plt.plot(ans.y[0],'r',label='R')
    plt.legend()

    plt.figure()
    plt.plot(ans.y[1],'r',label=r'$\frac{\partial R}{\partial r}$')

    plt.legend()



plt.show()