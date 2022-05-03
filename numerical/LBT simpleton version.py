from LTB_model_functions import *


def dSdt_dRdrdt(S,t,p):
    R, dRdr = S
    r, RR, EE, MM, dMMdr, dRRdr, dEEdr, G, rho_FLRW, r_b, n, m, A, H, Lamb = p

    ans_dRdr = ddRdrdt(r, RR, EE, MM, dMMdr, dRRdr, dEEdr, G, rho_FLRW, r_b, n, m, A, H, Lamb)
    ans_R = dRdt(r, RR, EE, MM, dMMdr, dRRdr, dEEdr, G, rho_FLRW, r_b, n, m, A, H, Lamb)
    
    return [ans_R, ans_dRdr]




# units for conversion
one_Mpc = 3.086e22 # m
one_Gy = 3.156e16 #s
one_solar_mass = 1.989e30 #kg

# values of the constants in the eq's
Lamb = 0
A = 1e-7
r_b = 5e10
n = 2
m = 2

H_SI = 68 # km/s/Mpc -> Mpc/Gyr
H = 68*1e3*one_Gy/one_Mpc # Mpc/Gy

G_SI = 6.67e-11  #  m^3/kg/s^2
G = 6.67e-11*one_solar_mass*one_Gy**2/(one_Mpc**3) # Mpc^3/M_o*Gy^2

rho_FLRW_SI =8.7e27 # kg/m^3 # This is the critical density Ryden in the back
rho_FLRW = 8.7e27*one_Mpc**3/one_solar_mass # M_o/Mpc^3

print(H,G,rho_FLRW)
a_i = 1/1100 #initial scale factor.

'''
Making a for loop on the radius, where we will integrate with the ODEint; on time;
for each loop
We have to have an initial R, which I will for now just set to zero
'''

dRRdr = a_i
MM = 1
EE = 0
dMMdr = 1
dRRdr = 1
dEEdr = 1
t_start= 1e3
t_end = 14e6
dt = 1e3
r= 0.1
RR = r*a_i
dr = 100

num_interations = 1
time_tot = np.linspace(t_start,t_end,round(t_end/10))

RR_results = np.zeros((num_interations,len(time_tot)))
RRR = a = [[0 for x in range(len(time_tot))] for x in range(m)]

for i in range(0,num_interations):
    ans = []
    #The initial conditions are found for each r, and used in the ODE int integration
    
    #RR += 1 # this matters a lot
    
    EE = E(r, RR,EE, MM, dMMdr, dRRdr, dEEdr, G, rho_FLRW, r_b, n, m, A, H, Lamb)
    #print(EE)
    MM = M(r, RR, EE, MM, dMMdr, dRRdr, dEEdr, G, rho_FLRW, r_b, n, m, A, H, Lamb)
    dEEdr = dEdr(r, RR, EE, MM, dMMdr, dRRdr, dEEdr, G, rho_FLRW, r_b, n, m, A, H, Lamb)
    dMMdr= dMdr(r, RR, EE, MM, dMMdr, dRRdr, dEEdr, G, rho_FLRW, r_b, n, m, A, H, Lamb)


    # The constants under integration
    args_list =[r, RR, EE, MM, dMMdr, dRRdr, dEEdr, G, rho_FLRW, r_b, n, m, A, H, Lamb]
    # the initial value for the function(s) that are being integrated
    init_cond_dRdt = [a_i*r, a_i]
    #making the time of integration
    #time_tot = np.linspace(t, 100,100)
    
    r += dr

    ans = scipy.integrate.odeint(dSdt_dRdrdt,
        t=time_tot,
        y0=init_cond_dRdt,
        args=(args_list,)
        )
    #RR_results[i] = ans[0] #Now we have our R


#print('RR_results =',RRR)
ans = ans.T
R_vec = ans[0]
print('len(R)=',len(R_vec))
dRdr_vec = ans[1]
print('len(dRdr)=',len(dRdr_vec))
print('R_max /R_min=',max(R_vec)/min(R_vec), '\n',
    'dRdr_max/dRdr_min =',max(dRdr_vec)/min(dRdr_vec))
#print('ans shape = ',ans)
#print(ans[0])
#print('time shape =',time_tot.shape)
plt.figure()
plt.plot(time_tot.T,R_vec,'-o',label='R(t,r)')
plt.plot(time_tot.T,dRdr_vec,'-o',label=r'$\frac{\partial R}{\partial r}$')
plt.legend()
plt.show()


