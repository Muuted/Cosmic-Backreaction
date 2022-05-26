from time import time
from LTB_model_functions import *
from Einstein_de_sitter_functions import *
from The_constants import *
import plotly as py
import plotly.express as px
import plotly.graph_objects as go

Lamb, A, r_b, n, m, H_0, H_i, G, rho_c0, rho_i0, a_i, t_i, t_0,c= func_constants()


dr = 1        # change of r
r_i = dr         # distance, we start at origo 
num_steps = 10000 # number of steps between t_start and t_end
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

Breakies = False
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



# Results for the Einstein de Sitter model 
a_ES, rho_EdS, time_vec = Einstein_de_sitter(num_of_steps=num_steps)

len_time = len(time_tot)

plt.figure()
rho_r = np.transpose(ans_rho)
plt.plot(radi_vec,rho_r[0]/rho_EdS[0],label=r'$\dfrac{\rho_{LTB}}{\rho_{EdS}}$(r,t $\approx$'+f'{round(time_tot[0]*1e4,2)}e-4)')
plt.plot(radi_vec,rho_r[int(len_time/2)]/rho_EdS[int(len_time/2)],label=r'$\dfrac{\rho_{LTB}}{\rho_{EdS}}$(r,t $\approx$'+f'{round(time_tot[int(len_time/2)],2)})')
plt.plot(radi_vec,rho_r[int(len_time)-1]/rho_EdS[int(len_time-1)],label=r'$\dfrac{\rho_{LTB}}{\rho_{EdS}}$(r,t $\approx$'+f'{round(time_tot[int(len_time-1)],2)})')
plt.xlabel('r [Mpc]')
plt.ylabel(r'$\frac{\rho_{LTB}}{\rho_{EdS}}$')
plt.title(r'$\dfrac{\rho(t_j,r)}{\rho_{EdS}(t_j)}$ at start, middle and end time, time in Gyr.'
         )
plt.legend()

plt.figure()
plt.plot(radi_vec,rho_r[0]/rho_EdS[0],label=r'$\dfrac{\rho_{LTB}}{\rho_{EdS}}$(r,t$\approx$'+f'{round(time_tot[0]*1e4,2)}e-4)')
plt.plot(radi_vec,rho_r[1]/rho_EdS[1],label=r'$\dfrac{\rho_{LTB}}{\rho_{EdS}}$(r,t$\approx$'+f'{round(time_tot[1]*1e3,2)}e-3)')
plt.plot(radi_vec,rho_r[2]/rho_EdS[2],label=r'$\dfrac{\rho_{LTB}}{\rho_{EdS}}$(r,t$\approx$'+f'{round(time_tot[2]*1e3,2)}e-3)')
plt.xlabel('r [Mpc]')
plt.ylabel(r'$\frac{\rho_{LTB}}{\rho_{EdS}}$')
plt.title(r'$\dfrac{\rho(t_j,r)}{\rho_{EdS}(t_j)}$ at 3 first time steps, time in Gyr')
plt.legend()



plt.figure()
rho_r = np.transpose(ans_rho)
ans_RR_trans = np.transpose(ans_RR)


RR_1 = ans_RR_trans[0]
RR_2 = ans_RR_trans[int(len_time/2)]
RR_3 = ans_RR_trans[int(len_time-1)]


plt.plot(RR_1,rho_r[0]/rho_EdS[0],label=r'$\dfrac{\rho_{LTB}}{\rho_{EdS}}$(r,t$\approx$'+f'{round(time_tot[0]*1e4,2)}e-4)')
plt.plot(RR_2,rho_r[int(len_time/2)]/rho_EdS[int(len_time/2)],label=r'$\dfrac{\rho_{LTB}}{\rho_{EdS}}$(r,t$\approx$'+f'{round(time_tot[int(len_time/2)],2)})')
plt.plot(RR_3,rho_r[int(len_time)-1]/rho_EdS[int(len_time-1)],label=r'$\dfrac{\rho_{LTB}}{\rho_{EdS}}$(r,t$\approx$'+f'{round(time_tot[int(len_time-1)],2)})')
plt.xlabel(r'$R(t_j,r)$ [Mpc]')
plt.ylabel(r'$\frac{\rho_{LTB}}{\rho_{EdS}}$')
plt.title(r'$\dfrac{\rho(t_j,r)}{\rho_{EdS}(t_j)}$ at start, middle and end time, time in Gyr.'
         )
plt.legend()


RR_1 = ans_RR_trans[0]
RR_2 = ans_RR_trans[1]
RR_3 = ans_RR_trans[2]
plt.figure()
plt.plot(RR_1,rho_r[0]/rho_EdS[0],label=r'$\dfrac{\rho_{LTB}}{\rho_{EdS}}$(r,t$\approx$'+f'{round(time_tot[0]*1e4,2)}e-4)')
plt.plot(RR_2,rho_r[1]/rho_EdS[1],label=r'$\dfrac{\rho_{LTB}}{\rho_{EdS}}$(r,t$\approx$'+f'{round(time_tot[1]*1e3,2)}e-3)')
plt.plot(RR_3,rho_r[2]/rho_EdS[2],label=r'$\dfrac{\rho_{LTB}}{\rho_{EdS}}$(r,t$\approx$'+f'{round(time_tot[2]*1e3,2)}e-3)')
plt.xlabel(r'$R(t_j,r)$ [Mpc]')
plt.ylabel(r'$\dfrac{\rho_{LTB}}{\rho_{EdS}}$')
plt.title(r'$\dfrac{\rho(t_j,r)}{\rho_{EdS}(t_j)}$ at 3 first time steps, time in Gyr')
plt.legend()
plt.show()

#------------------------------------------------------------ Finding volume element -----------------------
"""
Volume_ele = []

for j in range(0,len(time_tot)):
    V = 0
    j = 1
    for i in range(0,len(radi_vec)):
        if radi_vec[i] >= r_b:
            break
    
        V_dRdr_E = ans_dRdr[i][j]/(1+2*ans_E[i])
        V_R = ans_RR[i][j]

        V += 4*np.pi*V_dRdr_E*V_R**2*radi_vec[i]**2

    Volume_ele.append(V)


plt.figure()
plt.plot(time_tot,Volume_ele,label='Volume element')
plt.title('Volume element sum')
plt.xlabel('t [Gyr]')
plt.figure()
"""

#------------------------------------------------------------ Founding volume element -----------------------










plt.subplot(2,3,1)
for i in range(0,len(radi_vec),int(len(radi_vec)/4)):
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
#plt.show()







"""
plt.figure()
    plt.subplot(2,3,1)
    plt.plot(time_tot,ans_RR[len(ans_RR)-1]/radi_vec[len(radi_vec)-1], label=r'$\dfrac{R(t,r)}{r}$')
    plt.xlabel('t [Gyr]')
    plt.legend()


    plt.subplot(2,3,2)
    plt.plot(time_vec,a_ES,label=r'$a_{EdS}$')
    plt.xlabel('t [Gyr]')
    plt.legend()
    #plt.title('a(t) Ein de Sitter')

    plt.subplot(2,3,3)
    plt.plot(time_tot,ans_RR[len(ans_RR)-1]/radi_vec[len(radi_vec)-1], label=r'$\dfrac{R(t,r)}{r}$')
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
    plt.plot(time_vec,rho_EdS,label=r'$\rho_{EdS}$')

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



    rho_of_r1 = [[],[],[]]
    rho_of_r2 = []
    rho_of_r3 = []

    tmid =int(len(time_tot)/2)
    tend = len(time_tot)-1

    for i in range(0,len(radi_vec)):
        rho_of_r1[0].append(
            func_rho(radi_vec[i], ans_RR[i][0], ans_dMdr[i], ans_dRdr[i][0], rho_c0)/rho_EdS[0]
        )
        rho_of_r1[1].append(
            func_rho(radi_vec[i], ans_RR[i][tmid], ans_dMdr[i], ans_dRdr[i][tmid], rho_c0)/rho_EdS[tmid]
        )
        rho_of_r1[2].append(
            func_rho(radi_vec[i], ans_RR[i][tend], ans_dMdr[i], ans_dRdr[i][tend], rho_c0)/rho_EdS[tend]
        )

    print(np.shape(rho_of_r1))

    plt.figure()
    plt.plot(radi_vec,rho_of_r1[0],label=f'rho(r,t={time_tot[0]}')
    plt.plot(radi_vec,rho_of_r1[1],label=f'rho(r,t={time_tot[tmid]}')
    plt.plot(radi_vec,rho_of_r1[2],'--',label=f'rho(r,t={time_tot[tend]}')
    plt.legend()
    plt.title(r'$\rho(t_i,r)$/$\rho_{EdS}[0]$')
    plt.xlabel('r [Mpc]')


    plt.figure()
    for i in range(0,len(radi_vec),10):
        plt.plot(time_tot,
        func_rho(radi_vec[i], ans_RR[i], ans_dMdr[i], ans_dRdr[i], rho_c0)/rho_EdS
        ,'--',label=f'rho(t,r={radi_vec[i]})/rho_EdS'
        )

    #.xlim(0.5,time_tot[len(time_tot)-1])
    #plt.ylim(1-0.00015, 1.00015)
    plt.title(r'evolution of $\dfrac{\rho(t,r_i)}{\rho_{EdS}}$ ylim(1-0.00015, 1.00015)')
    plt.legend()

    plt.figure()
    for i in range(0,len(radi_vec),10):
        plt.plot(time_tot,
        func_rho(radi_vec[i], ans_RR[i], ans_dMdr[i], ans_dRdr[i], rho_c0)/rho_EdS
        ,'--',label=f'rho(t,r={radi_vec[i]})/rho_EdS'
        )

    #plt.ylim(0.99995, 1.00005)
    plt.title(r'evolution of $\dfrac{\rho(t,r_i)}{\rho_{EdS}}$ ylim = non')
    plt.legend()


    plt.figure()
    for i in range(0,len(radi_vec),10):
        plt.plot(time_tot,
        ans_RR[i]/radi_vec[i]
        ,'-',label=f'R(t,r={radi_vec[i]})'
        )
    plt.plot(time_vec,a_ES,'--',label=f'$a_(EdS)')
    plt.title('Evolution of R/r at different r')
    plt.xlabel('t [Gyr]')
    plt.legend()


    plt.figure()
    for i in range(0,len(radi_vec),10):
        plt.plot(time_tot,
        ans_dRdr[i]
        ,'-',label=f'dRdr(t,r={radi_vec[i]})'
        )
    plt.title('Evolution of dRdr at different r')
    plt.legend()


    avg_R = []
    avgR = 0
    for j in range(len(ans_RR[0])):
        for i in range(len(radi_vec)):
            avgR += ans_RR[i][j]/(len(radi_vec)*radi_vec[i])
        avg_R.append(avgR)
    print(len(avg_R),len(ans_RR[0]))

    plt.figure()
    plt.plot(time_tot,avg_R,label=r'$R_{avg}$')
    plt.title(r'Averaging all $\dfrac{R(t,r_i)}{r_i}$ gives = $R_{avg}$')
    plt.ylabel('t [Gyr]')
    plt.legend()


 
    y = time_tot
    x = radi_vec
    z = ans_rho

    #fig = go.Figure(data=[go.Surface(z=z, x=x, y=y)])
    #fig.update_layout(title='Mt Bruno Elevation', autosize=False,
                    #   width=500, height=500,
                    #  margin=dict(l=65, r=50, b=65, t=90))
    #fig.show()

  

    plot_rho = [[],[],[],[],[]]
    r_pos_vec = []
    k = 0
    for i in range(0,len(radi_vec),10):
        if k >5 :
            break
        r_pos_vec.append(radi_vec[i])
        for j in range(0,len(ans_RR[0])-1):
            plot_rho[k].append(
                func_rho(radi_vec[i], ans_RR[i][j], ans_dMdr[i], ans_dRdr[i][j], rho_c0)/rho_EdS[j]
                )
        k+=1
        if k >5 :
            break


    plot1 = plot_rho[0]
    plot2 = plot_rho[1]
    plot3 = plot_rho[2]
    plot4 = plot_rho[3]
    plot5 = plot_rho[4]

    fig2 = go.Figure()
    fig2.add_trace(
        go.Scatter(x=time_tot,y=plot1
        ,mode='lines',name=f'rho(t,r={r_pos_vec[0]})'
        )
        )
    fig2.add_trace(
        go.Scatter(x=time_tot,y=plot2
        ,mode='lines',name=f'rho(t,r={r_pos_vec[1]})'
        )
    )
    fig2.add_trace(
    go.Scatter(x=time_tot,y=plot3
    ,mode='lines',name=f'rho(t,r={r_pos_vec[2]})'
    )
    )
    fig2.add_trace(
    go.Scatter(x=time_tot,y=plot4
    ,mode='lines+markers',name=f'rho(t,r={r_pos_vec[3]})'
    )
    )
    fig2.add_trace(
    go.Scatter(x=time_tot,y=plot5
        ,name=f'rho(t,r={r_pos_vec[4]})',line=dict(dash='dash',color='black')
    )
    )
    fig2.update_layout(
        yaxis_range=[0.99994,1.00004],
        title=r'$\rho(t,r_i) \text{  time evolution for different starting positions}$'
        )
    fig2.show()
   

    rhos = []
    print(len(radi_vec))
    for i in range(0,len(radi_vec)):
        rhos.append(
            func_rho(radi_vec[i],a_i*radi_vec[i],ans_dMdr[i],a_i,rho_c0)/func_rho_Ein_Sitter(a_i)
        )
    print(len(rhos))
    plt.figure()
    plt.plot(radi_vec,rhos)

    plt.show()

    r = 2
    RR = r*a_i
    dRRdr = a_i
    EE = func_E(r, r_b, n, m, A)
    dEEdr = func_dEdr(r, r_b, n, m, A)
    MM = func_M(r, EE, G, rho_c0, a_i, H_i, c)
    dMMdr = func_dMMdr(r, EE, dEEdr, rho_c0, H_i, a_i, c) 

    print('M=',MM/1e12)
    print('dMdr=',dMMdr/1e12)


    print('RR=',ans_RR[3][0],RR,'\n','r=',radi_vec[3]/1100,'\n','t=',time_tot[len(ans_RR[0])-1])
    print('dRdr=',ans_dRdr[3][len(ans_RR[0])-1],'\n','r=',radi_vec[3],'\n','t=',time_tot[len(ans_RR[0])-1])


"""