from LTB_model_functions import *


def dSdt_dRdrdt(S,t,p):
    dRdr = S
    r,dMMdr, dRRdr, dEEdr, G, rho_i, r_b, n, m, A, H, Lamb = p

    ans_dRdrdt = ddRdrdt(r, RR, EE, MM, dMMdr, dRRdr, dEEdr, G, rho_i, r_b, n, m, A, H, Lamb)

    return ans_dRdrdt

def dSdt_dRdt(S,t,p):
    R = S
    r, RR, EE, MM, dMMdr, dRRdr, dEEdr, G, rho_i, r_b, n, m, A, H, Lamb = p

    ans_dRdt = dRdt(r, RR, EE, MM, dMMdr, dRRdr, dEEdr, G, rho_i, r_b, n, m, A, H, Lamb)

    return ans_dRdt


# values of the constants in the eq's
Lamb = 0
A = 1e-7
r_b = 200 #5e10
n = 2
m = 2
H = 1
G = 6.67e-11  # Nm**2/kg**2 = m**3/kg*s**2
rho_i = 1  # 1e-26 #kg/m^3
args_list = (G, rho_i, r_b, n, m, A, H, Lamb)

'''
Making a for loop on the radius, where we will integrate with the ODEint; on time;
for each loop
We have to have an initial R, which I will for now just set to zero
'''
RR = 1
dRRdr = 1
MM = 1
EE = 0
dMMdr = 100
dRRdr = 1
dEEdr = 1
t= 0
r= r_b
dr = 100

num_interations = 1
for i in range(0,num_interations):
    
    #The initial conditions are found for each r, and used in the ODE int integration
    
    #RR += 1 # this matters a lot
    EE = E(r, RR,EE, MM, dMMdr, dRRdr, dEEdr, G, rho_i, r_b, n, m, A, H, Lamb)
    #print(EE)
    MM = M(r, RR, EE, MM, dMMdr, dRRdr, dEEdr, G, rho_i, r_b, n, m, A, H, Lamb)
    dEEdr = dEdr(r, RR, EE, MM, dMMdr, dRRdr, dEEdr, G, rho_i, r_b, n, m, A, H, Lamb)
    dMMdr= dMdr(r, RR, EE, MM, dMMdr, dRRdr, dEEdr, G, rho_i, r_b, n, m, A, H, Lamb)


    # The constants under integration
    args_list =[r, RR, EE, MM, dMMdr, dRRdr, dEEdr, G, rho_i, r_b, n, m, A, H, Lamb]
    # the initial value for the function(s) that are being integrated
    init_cond_dRdt = [0]
    #making the time of integration
    time_tot = np.linspace(t, 100,100)
    
    r += dr

    ans = scipy.integrate.odeint(dSdt_dRdt,
        t=time_tot,
        y0=init_cond_dRdt,
        args=(args_list,)
        )
    RR = ans[0].T #Now we have our R

    '''ans2 = scipy.integrate.odeint(dSdt_dRdrdt,
        t=time_tot,
        args=args_list,
        )'''

ans = ans.T
print('ans shape = ',ans[0].shape)
#print(ans[0])
print('time shape =',time_tot.shape)
plt.figure()
plt.plot(time_tot.T,ans[0])
plt.show()

