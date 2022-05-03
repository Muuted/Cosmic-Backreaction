from LTB_model_functions import *


def dSdt_dRdrdt(S,t,p):
    R, dRdr = S
    r, RR, EE, MM, dMMdr, dRRdr, dEEdr, G, rho_i, r_b, n, m, A, H, Lamb = p

    ans_dRdr = ddRdrdt(r, RR, EE, MM, dMMdr, dRRdr, dEEdr, G, rho_i, r_b, n, m, A, H, Lamb)
    ans_R = dRdt(r, RR, EE, MM, dMMdr, dRRdr, dEEdr, G, rho_i, r_b, n, m, A, H, Lamb)
    
    return [ans_dRdr, ans_R]

def dSdt_dRdt(S,t,p):
    R= S
    r, RR, EE, MM, dMMdr, dRRdr, dEEdr, G, rho_i, r_b, n, m, A, H, Lamb = p

    ans_dRdt = dRdt(r, RR, EE, MM, dMMdr, dRRdr, dEEdr, G, rho_i, r_b, n, m, A, H, Lamb)

    return ans_dRdt


# values of the constants in the eq's
Lamb = 0
A = 1e-7
r_b = 5e10
n = 2
m = 2
H = 70e3 # m/s/Mpc
G = 6.67e-11  # Nm**2/kg**2 = m**3/kg*s**2
rho_i =  1e-26 #kg/m^3
a_i = 1/1100
args_list = (G, rho_i, r_b, n, m, A, H, Lamb)

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
t= 0
r= 0.1
RR = r*a_i
dr = 100

num_interations = 1
time_tot = np.linspace(0, 100000,100)
RR_results = np.zeros((num_interations,len(time_tot)))
RRR = a = [[0 for x in range(len(time_tot))] for x in range(m)]

for i in range(0,num_interations):
    ans = []
    #The initial conditions are found for each r, and used in the ODE int integration
    
    #RR += 1 # this matters a lot
    rho_i = rho(r, RR,EE, MM, dMMdr, dRRdr, dEEdr, G, rho_i, r_b, n, m, A, H, Lamb)
    EE = E(r, RR,EE, MM, dMMdr, dRRdr, dEEdr, G, rho_i, r_b, n, m, A, H, Lamb)
    #print(EE)
    MM = M(r, RR, EE, MM, dMMdr, dRRdr, dEEdr, G, rho_i, r_b, n, m, A, H, Lamb)
    dEEdr = dEdr(r, RR, EE, MM, dMMdr, dRRdr, dEEdr, G, rho_i, r_b, n, m, A, H, Lamb)
    dMMdr= dMdr(r, RR, EE, MM, dMMdr, dRRdr, dEEdr, G, rho_i, r_b, n, m, A, H, Lamb)


    # The constants under integration
    args_list =[r, RR, EE, MM, dMMdr, dRRdr, dEEdr, G, rho_i, r_b, n, m, A, H, Lamb]
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
print('ans shape = ',ans[0].shape)
#print(ans[0])
print('time shape =',time_tot.shape)
plt.figure()
plt.plot(time_tot.T,ans[0],'-o')
plt.show()

