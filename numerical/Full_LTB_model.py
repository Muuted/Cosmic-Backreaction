from logging import PlaceHolder
from LTB_model_functions import *
from Einstein_de_sitter_functions import *
from The_constants import *
import plotly as py
import plotly.graph_objects as go

Lamb, A, r_b, n, m, H_0, H_i, G, rho_c0, rho_i0, a_i, t_i, t_0,c= func_constants()

dt = 1e3        # time step
dr = 1        # change of r
r_i = dr         # distance, we start at origo 
num_steps = 10000 # number of steps between t_start and t_end
num_iterations = int(r_b/dr) #number of r's

# Our time vector for the integration
time_tot = np.linspace(t_i,t_0,num_steps)
radi_vec = np.linspace(r_i,r_b,num_iterations)

# Creating the lists for the data at different r values
ans_RR = [[] for x in range(num_iterations)]
ans_dRdr = [[] for x in range(num_iterations)]
ans_rho = [[] for x in range(num_iterations)]

ans_M = []
ans_dMdr = []
ans_dEdr = []
ans_E = []

Breakies = False
for i in range(0,num_iterations): # This loop makes it so that we iterate over r
    # increasing r for each iteration
    r = radi_vec[i]
    
    # Initial condition for the functions
    RR = r*a_i
    dRRdr = a_i

    EE = func_E(r, r_b, n, m, A)
    dEEdr = func_dEdr(r, r_b, n, m, A)
    MM = func_M(r, EE, G, rho_c0, a_i, H_i, c)
    dMMdr = func_dMMdr(r, EE, dEEdr, rho_c0, H_i, a_i, c) 


    ans_E.append(EE)
    ans_dEdr.append(dEEdr)
    ans_M.append(MM)
    ans_dMdr.append(dMMdr)

    if i ==0:
        print('E=',EE,'\n','dEdr=',dEEdr,'\n','M=',MM,'\n','dMdr=',dMMdr,'\n','r=',r)
    dRdt = func_dRdt(r, RR, EE, MM, G, c)
    dRdrdt = func_dRdrdt(r, RR, EE, MM, dMMdr, dRRdr, dEEdr, G, c)
    #print('inner loop',dRdt,dRdrdt)

    #The initial conditions are found for each r, and used in the ODE int integration
    init_cond_dRdt = [RR, dRRdr]
    args_list =[r, EE, dEEdr, MM,dMMdr, G,c]

    #print(init_cond_dRdt)
    ans_odeint = scipy.integrate.odeint(func_LTB_dSdt, 
            y0=init_cond_dRdt, 
            t=time_tot,
            args=(args_list,)
            )
    # The results from our odeint  above
    ans_odeint = ans_odeint.T
    
    ans_RR[i] = ans_odeint[0]
    ans_dRdr[i] = ans_odeint[1]

    for j in range(0,len(ans_RR[i])-1):
        rho_valu = func_rho(r,ans_RR[i][j],dMMdr,ans_dRdr[i][j],rho_c0)
        ans_rho[i].append(rho_valu)
        
        if str(rho_valu) == 'nan' or str(ans_RR[i][j])=='nan':
            print('RR = ',ans_RR[i][j],'\n'
            ,'dRdr =', ans_dRdr[i][j], '\n'
            ,f'E(r={r})=',EE,'\n'
            ,f'dE(r={r})dr=',dEEdr,'\n'
            ,f'M(r={r})=',MM,'\n'
            ,f'dM(r={r})dr=',dMMdr,'\n'
            ,f'rho(r={r})dr=',ans_rho[i][j]
            )
            Breakies = True
            break
    if Breakies == True:
        break


"""
# Results for the Einstein de Sitter model 
a_ES, rho, rho_ES, time_vec = Einstein_de_sitter(num_of_steps=num_steps)
ans_a_ES = rho_ES

plt.figure()
plt.subplot(2,3,1)
plt.plot(time_tot,ans_RR[0]/r, label=r'$\dfrac{R(t,r)}{r}$')
plt.xlabel('t [Gyr]')
plt.legend()


plt.subplot(2,3,2)
plt.plot(time_vec,a_ES,label=r'$a_{EdS}$')
plt.xlabel('t [Gyr]')
plt.legend()
#plt.title('a(t) Ein de Sitter')

plt.subplot(2,3,3)
plt.plot(time_tot,ans_RR[0]/r,label=r'$\dfrac{R(t,r)}{r}$')
plt.plot(time_vec,a_ES,'--',label=r'$a_{EdS}$')
plt.xlabel('t [Gyr]')
#plt.title('R vs a EdS')
plt.legend()

plt.subplot(2,3,4)
plt.plot(time_tot,ans_dRdr[0],label=r'$\dfrac{\partial R(t,r)}{\partial r}$')
plt.xlabel('t [Gyr]')
plt.legend()
#plt.title('dRdr')


plt.subplot(2,3,5)
plt.plot(time_tot,func_rho(r, ans_RR[0], dMMdr, ans_dRdr[0], rho_c0),'-o',label=r'$\rho_{LTB}$(t,r)')
plt.plot(time_vec,rho,label=r'$\rho_{EdS}$')

plt.xlabel('t [Gyr]')
plt.legend()
#plt.title('rho')

E_vec = []
for i in range(0,len(radi_vec)):
    E_vec.append(func_E(radi_vec[i],r_b,n,m,A))

plt.subplot(2,3,6)
plt.plot(radi_vec,E_vec,label=r'E(r)')
plt.xlabel('r [Mpc]')
plt.legend()

"""

print(np.shape(ans_RR))
print(np.shape(np.transpose(ans_RR)))
print(np.shape(radi_vec))
plt.figure()
plt.subplot(2,3,1)
plt.plot(radi_vec,ans_E,label='E(r)')
plt.xlabel('r [Mpc]')
plt.legend()

plt.subplot(2,3,2)
plt.plot(radi_vec,ans_M,label='M(r)')
plt.xlabel('r [Mpc]')
plt.legend()

plt_RR =np.transpose(ans_RR)
plt.subplot(2,3,3)
plt.plot(radi_vec,plt_RR[0],label='R(r)')
plt.xlabel('r [Mpc]')
plt.legend()

plt.subplot(2,3,4)
plt.plot(radi_vec,ans_dEdr,label='dEdr')
plt.xlabel('r [Mpc]')
plt.legend()


plt.subplot(2,3,5)
plt.plot(radi_vec,ans_dMdr,label='dMdr')
plt.xlabel('r [Mpc]')
plt.legend()

plt_dRdr = np.transpose(ans_dRdr)
plt.subplot(2,3,6)
plt.plot(radi_vec,plt_dRdr[0],label='dRdr(r)')
plt.xlabel('r [Mpc]')
plt.legend()



plt.figure()
i = 0
plt.subplot(2,3,1)
plt.plot(time_tot,func_rho(radi_vec[i], ans_RR[i], ans_dMdr[i], ans_dRdr[i], rho_c0),'-o',label=f'rho(t,r={radi_vec[i]})')
plt.xlabel('t [Gyr]')

plt.legend()

i = 10
plt.subplot(2,3,2)
plt.plot(time_tot,func_rho(radi_vec[i], ans_RR[i], ans_dMdr[i], ans_dRdr[i], rho_c0),'-o',label=f'rho(t,r={radi_vec[i]})')
plt.xlabel('t [Gyr]')
plt.legend()

i = 20
plt.subplot(2,3,3)
plt.plot(time_tot,func_rho(radi_vec[i], ans_RR[i], ans_dMdr[i], ans_dRdr[i], rho_c0),'-o',label=f'rho(t,r={radi_vec[i]})')
plt.xlabel('t [Gyr]')
plt.legend()

i = 40
plt.subplot(2,3,4)
plt.plot(time_tot,func_rho(radi_vec[i], ans_RR[i], ans_dMdr[i], ans_dRdr[i], rho_c0),'-o',label=f'rho(t,r={radi_vec[i]})')
plt.xlabel('t [Gyr]')
plt.legend()

i = 49
plt.subplot(2,3,5)
plt.plot(time_tot,func_rho(radi_vec[i], ans_RR[i], ans_dMdr[i], ans_dRdr[i], rho_c0),'-o',label=f'rho(t,r={radi_vec[i]})')
plt.xlabel('t [Gyr]')
plt.legend()

plt.subplot(2,3,6)
#plt.plot(radi_vec,func_rho(radi_vec, ans_RR[0][0], ans_dMdr[0], ans_dRdr[0][0], rho_c0),'-o',label=f'rho(t=t_i,r)')
#plt.plot(radi_vec,ans_dMdr)
plt.xlabel('r [Mpc]')
plt.legend()



plt.figure()
for i in range(0,len(radi_vec),10):
    plt.plot(time_tot,func_rho(radi_vec[i], ans_RR[i], ans_dMdr[i], ans_dRdr[i], rho_c0),'-o',label=f'rho(t,r={radi_vec[i]})')

#plt.xlim(0.0011870,time_tot[len(time_tot)-1])
#plt.ylim(-2,1e14)
plt.legend()

plt.show()

y = time_tot
x = radi_vec
z = ans_rho

fig = go.Figure(data=[go.Surface(z=z, x=x, y=y)])
fig.update_layout(title='Mt Bruno Elevation', autosize=False,
                    width=500, height=500,
                   margin=dict(l=65, r=50, b=65, t=90))
#fig.show()

