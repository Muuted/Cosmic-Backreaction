# ------------------------ Full_LTB_model.py ----------------------------
#from time import time
from matplotlib.markers import MarkerStyle
from LTB_model_functions import *
from Einstein_de_sitter_functions import *
from The_constants import *
import plotly as py
import plotly.express as px
import plotly.graph_objects as go

Lamb, A, r_b, n, m, H_0, H_i, G, rho_c0, rho_i0, a_i, t_i, t_0,c= func_constants()

dr = 1        # change of r
r_i = dr         # distance, we start at origo 
num_steps = 1000 # number of steps between t_start and t_end
num_iterations = 100#int(r_b/dr) #number of r's

# Our time vector for the integration
time_tot = np.linspace(t_i,t_0,num_steps)
radi_vec = np.linspace(r_i,r_b*1.2,num_iterations)


# Creating the lists for the data at different r values
ans_RR = [[] for x in range(num_iterations)]
ans_dRdr = [[] for x in range(num_iterations)]
ans_rho = [[] for x in range(num_iterations)]

ans_M = []
ans_dMdr = []
ans_dEdr = []
ans_E = []

for i in range(0,num_iterations): # This loop makes it so that we iterate over r
    # Getting the r for each iteration.
    r = radi_vec[i]

    # Calculating the initial conditions
    RR = r*a_i
    dRRdr = a_i
    EE = func_E(r, r_b, n, m, A)
    dEEdr = func_dEdr(r, r_b, n, m, A)
    MM = func_M(r, EE, G, rho_c0, a_i, H_i, c)
    dMMdr = func_dMMdr(r, EE, dEEdr, rho_c0, H_i, a_i, c) 
    
    #The initial conditions are found for each r, and used in the ODE int integration
    init_cond_dRdt = [RR, dRRdr]
    args_list =[EE, dEEdr, MM,dMMdr, G,c]

    ans_odeint = scipy.integrate.odeint(func_LTB_dSdt, 
            y0=init_cond_dRdt, 
            t=time_tot,
            args=(args_list,)
            )
    # The results from our odeint  above

    ans_odeint = ans_odeint.T
    
    ans_RR[i] = ans_odeint[0]
    ans_dRdr[i] = ans_odeint[1]

    
    ans_E.append(EE)
    ans_dEdr.append(dEEdr)
    ans_M.append(MM)
    ans_dMdr.append(dMMdr)

    for j in range(0,len(ans_RR[i])):
        rho_valu = func_rho(ans_RR[i][j],dMMdr,ans_dRdr[i][j],rho_c0)
        ans_rho[i].append(rho_valu)
        


# Results for the Einstein de Sitter model 
a_ES, rho_EdS, time_vec = Einstein_de_sitter(num_of_steps=num_steps)


#------------------------------------------------------------ Finding volume element -----------------------

Volume_LTB = []
Volume_EdS = []
for j in range(0,len(time_tot)):
    V_EdS = 0
    V_LTB = 0
    V_LTB_beackreac = 0
    for i in range(0,len(radi_vec)):
    
        V_dRdr_E = ans_dRdr[i][j]/np.sqrt(1+2*ans_E[i])
        V_R = ans_RR[i][j]
        
        V_LTB += 4*np.pi*V_dRdr_E*V_R**2*(max(radi_vec)/(len(radi_vec)))

    V_EdS = (4*np.pi/3)*(a_ES[j]*max(radi_vec))**3

    Volume_LTB.append(V_LTB)
    Volume_EdS.append(V_EdS)
   

plt.figure()
plt.plot(time_tot,Volume_LTB,label=r'$V_{LTB,\mathcal{D}}$')
plt.plot(time_tot,Volume_EdS,label=r'$V_{EdS,\mathcal{D}}$')
plt.title(r'Volume fraction $V_{\mathcal{D}}$ from r = 0 to '+f'r={max(radi_vec)} where '+r'$r_b$='+f'{r_b}',fontsize = 15)
plt.xlabel('t [Gyr]',fontsize=15)
plt.ylabel(r'V [$Mpc^3$]',fontsize = 15)
plt.legend(fontsize = 15)

#------------------------------------------------------------ Founding volume element -----------------------

k = 0
pos_vecs = []
i_pos = []

for i in range(0,len(radi_vec)):
    if k < 1 and radi_vec[i] < 10:
        pos_vecs.append(radi_vec[i])
        i_pos.append(i)
        k += 1
    if k < 2 and 20 < radi_vec[i] < 25:
        pos_vecs.append(radi_vec[i])
        i_pos.append(i)
        k += 1
    if k < 3 and 40 < radi_vec[i] < 42:
        pos_vecs.append(radi_vec[i])
        i_pos.append(i)
        k += 1
    if k < 4 and r_b*1.1 < radi_vec[i]:
        pos_vecs.append(radi_vec[i])
        i_pos.append(i)
        k += 1



plt.figure()
plt.plot(time_tot,ans_RR[i_pos[0]]/radi_vec[i_pos[0]],marker='.',markersize=5,label=r'R(t,$r_1$ $\approx$ '+f'{int(pos_vecs[0])})'+r'/$~r_1$')
plt.plot(time_tot,ans_RR[i_pos[1]]/radi_vec[i_pos[1]],marker='.',markersize=5,label=r'R(t,$r_2$ $\approx$ '+f'{int(pos_vecs[1])})'+r'/$~r_2$')
plt.plot(time_tot,ans_RR[i_pos[2]]/radi_vec[i_pos[2]],marker='.',markersize=5,label=r'R(t,$r_3$ $\approx$ '+f'{int(pos_vecs[2])})'+r'/$~r_3$')
plt.plot(time_tot,ans_RR[i_pos[3]]/radi_vec[i_pos[3]],marker='.',markersize=5,label=r'R(t,$r_4$ $\approx$ '+f'{int(pos_vecs[3])})'+r'/$~r_4$')
plt.plot(time_tot,a_ES,'--k',label=r'$a_{EdS}(t)$')
plt.title(r'Proper distance $R(t,r_j)$ at different radial positions',fontsize=15)
plt.xlabel('t [Gyr]',fontsize=20)
plt.ylabel(' R(t,r) [Mpc]',fontsize=20)
plt.legend(fontsize=15)



len_time = len(time_tot)

plt.figure()
rho_r = np.transpose(ans_rho)
plt.plot(radi_vec,rho_r[0]/rho_EdS[0],label=r'$\dfrac{\rho_{LTB}}{\rho_{EdS}}$(r,t $\approx$'+f'{round(time_tot[0]*1e4,2)}e-4)')
plt.plot(radi_vec,rho_r[int(len_time/2)]/rho_EdS[int(len_time/2)],label=r'$\dfrac{\rho_{LTB}}{\rho_{EdS}}$(r,t $\approx$'+f'{round(time_tot[int(len_time/2)],2)})')
plt.plot(radi_vec,rho_r[int(len_time)-1]/rho_EdS[int(len_time-1)],label=r'$\dfrac{\rho_{LTB}}{\rho_{EdS}}$(r,t $\approx$'+f'{round(time_tot[int(len_time-1)],2)})')
plt.xlabel('r [Mpc]',fontsize=20)
plt.ylabel(r'$\frac{\rho_{LTB}}{\rho_{EdS}}$',fontsize=20)
plt.title(r'$\dfrac{\rho(t_j,r)}{\rho_{EdS}(t_j)}$ at start, middle and end time, time in Gyr.',fontsize=15)
plt.legend(fontsize=15,loc='lower right')

plt.figure()
plt.plot(radi_vec,rho_r[0]/rho_EdS[0],label=r'$\dfrac{\rho_{LTB}}{\rho_{EdS}}$(r,t$\approx$'+f'{round(time_tot[0]*1e4,2)}e-4)')
plt.plot(radi_vec,rho_r[1]/rho_EdS[1],label=r'$\dfrac{\rho_{LTB}}{\rho_{EdS}}$(r,t$\approx$'+f'{round(time_tot[1]*1e3,2)}e-3)')
plt.plot(radi_vec,rho_r[2]/rho_EdS[2],label=r'$\dfrac{\rho_{LTB}}{\rho_{EdS}}$(r,t$\approx$'+f'{round(time_tot[2]*1e3,2)}e-3)')
plt.xlabel('r [Mpc]',fontsize=20)
plt.ylabel(r'$\frac{\rho_{LTB}}{\rho_{EdS}}$',fontsize=20)
plt.title(r'$\dfrac{\rho(t_j,r)}{\rho_{EdS}(t_j)}$ at 3 first time steps, time in Gyr',fontsize=15)
plt.legend(fontsize=15,loc='lower right')



plt.figure()
rho_r = np.transpose(ans_rho)
ans_RR_trans = np.transpose(ans_RR)


RR_1 = ans_RR_trans[0]
RR_2 = ans_RR_trans[int(len_time/2)]
RR_3 = ans_RR_trans[int(len_time-1)]


plt.plot(RR_1,rho_r[0]/rho_EdS[0],marker='.',markersize=7,label=r'$\dfrac{\rho_{LTB}}{\rho_{EdS}}$(r,t$\approx$'+f'{round(time_tot[0]*1e4,2)}e-4)')
plt.plot(RR_2,rho_r[int(len_time/2)]/rho_EdS[int(len_time/2)],marker='.',markersize=7,label=r'$\dfrac{\rho_{LTB}}{\rho_{EdS}}$(r,t$\approx$'+f'{round(time_tot[int(len_time/2)],2)})')
plt.plot(RR_3,rho_r[int(len_time)-1]/rho_EdS[int(len_time-1)],marker='.',markersize=7,label=r'$\dfrac{\rho_{LTB}}{\rho_{EdS}}$(r,t$\approx$'+f'{round(time_tot[int(len_time-1)],2)})')
plt.xlabel(r'$R(t_j,r)$ [Mpc]',fontsize=20)
plt.ylabel(r'$\dfrac{\rho_{LTB}}{\rho_{EdS}}$',fontsize=20)
plt.title(r'$\dfrac{\rho(t_j,r)}{\rho_{EdS}(t_j)}$ at start, middle and end time, time in Gyr.', fontsize=15
         )
plt.legend(fontsize=15,title_fontsize=1,loc='upper left')



RR_1 = ans_RR_trans[0]
RR_2 = ans_RR_trans[1]
RR_3 = ans_RR_trans[2]
plt.figure()
plt.plot(RR_1,rho_r[0]/rho_EdS[0],marker='.',markersize=5,label=r'$\dfrac{\rho_{LTB}}{\rho_{EdS}}$(r,t$\approx$'+f'{round(time_tot[0]*1e4,2)}e-4)')
plt.plot(RR_2,rho_r[1]/rho_EdS[1],marker='.',markersize=5,label=r'$\dfrac{\rho_{LTB}}{\rho_{EdS}}$(r,t$\approx$'+f'{round(time_tot[1]*1e3,2)}e-3)')
plt.plot(RR_3,rho_r[2]/rho_EdS[2],marker='.',markersize=5,label=r'$\dfrac{\rho_{LTB}}{\rho_{EdS}}$(r,t$\approx$'+f'{round(time_tot[2]*1e3,2)}e-3)')
plt.xlabel(r'$R(t_j,r)$ [Mpc]',fontsize=20)
plt.ylabel(r'$\dfrac{\rho_{LTB}}{\rho_{EdS}}$',fontsize=20)
plt.title(r'$\dfrac{\rho(t_j,r)}{\rho_{EdS}(t_j)}$ at 3 first time steps, time in Gyr',fontsize=15)
plt.legend(fontsize= 15, loc ='lower right')



plt.subplot(2,3,1)
for i in range(0,len(radi_vec),int(len(radi_vec)/5)):
    plt.plot(time_tot,
    ans_RR[i]/radi_vec[i]
    ,'-',label=f'R(t,r~{int(radi_vec[i])})'
    )
#plt.plot(time_tot,ans_RR[len(radi_vec[0])-1]/radi_vec[len(radi_vec[0])-1],'-',label=f'R(t,r={radi_vec[len(radi_vec[0])-1]})')
plt.plot(time_vec,a_FLRW_lim(time_tot,t_0),'--',label=r'$a_{EdS}$')
plt.title('Evolution of R/r at different r')
plt.xlabel('t [Gyr]')
plt.legend()


plt.subplot(2,3,2)
t = 0
plt.plot(time_tot,ans_rho[0]/rho_EdS[0])
plt.plot(time_tot,ans_rho[25]/rho_EdS[1])
plt.plot(time_tot,ans_rho[49]/rho_EdS[2])
plt.xlabel('t [Gyr]')



plt.subplot(2,3,3)
rho_r = np.transpose(ans_rho)
plt.plot(radi_vec,rho_r[0]/rho_EdS[0],label='rho(r,t_i)')
plt.plot(radi_vec,rho_r[1]/rho_EdS[1],label='rho(r,t_mid)')
plt.plot(radi_vec,rho_r[2]/rho_EdS[2],label='rho(r,t_0)')
plt.xlabel('r [Mpc]')
plt.legend()



plt.subplot(2,3,4)
#plt.plot(time_vec,rho_EdS,label='rho_EdS')
#plt.plot(radi_vec,ans_dRdr,label='dRdr')
plt.plot(radi_vec,ans_RR)

plt.legend()
plt.subplot(2,3,5)
plt.plot(radi_vec,ans_M,label='M')
plt.legend()
plt.xlabel('r [Mpc]')

plt.subplot(2,3,6)
plt.plot(radi_vec,ans_dMdr,label='dMdr')
plt.legend()
plt.xlabel('r [Mpc]')

#plt.show()


fig1 = go.Figure()
fig1.add_trace(
    go.Scatter(x=radi_vec,y=rho_r[0]/rho_EdS[0]
    ,mode='lines',name=f'rho(t=~{int(time_tot[0])},r)'
        )
    )
fig1.update_layout(
    title=r'$\rho(t_i,r)/\rho_{EdS} $s'
    , xaxis_title='r [Mpc]'
)

fig2 = go.Figure()
fig2.add_trace(
    go.Scatter(x=radi_vec,y=rho_r[0]/rho_EdS[0]
    ,mode='lines',name=f'rho(t=~{int(time_tot[0])},r)'
        )
    )

fig2.add_trace(
    go.Scatter(x=radi_vec,y=rho_r[1]/rho_EdS[1]#/max(rho_r[1])
    ,mode='lines',name=f'rho(t=~{int(time_tot[1])},r)'
        )   
    )
fig2.add_trace(
    go.Scatter(x=radi_vec,y=rho_r[2]/rho_EdS[2]#/max(rho_r[2])
    ,mode='lines',name=f'rho(t=~{int(time_tot[2])},r)'
        )
    )

fig2.update_layout(
    title=r'$\rho(t,r)/\rho_{EdS} $s'
    , xaxis_title='r [Mpc]'
)
#fig2.show()
#fig1.show()

print(np.shape(ans_RR))
delta_R = []
delta_R.append(0)
for i in range(len(radi_vec)-1):
    delta_R.append(
        ans_RR[i+1][len(time_tot)-1] - ans_RR[i][len(time_tot)-1]
    )

plt.figure()
plt.plot(radi_vec,delta_R)


plt.show()




