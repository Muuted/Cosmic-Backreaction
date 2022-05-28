import numpy as np
import matplotlib.pyplot as plt
import scipy as sp


def v(theta, f_2):
    numerator = f_2 * np.pi ** 3 * (1 - np.cos(theta)) ** 3
    denominator = 8 * (1 - f_2) * (theta - np.sin(theta)) ** 3 + f_2 * np.pi ** 3 * (1 - np.cos(theta)) ** 3
    return numerator / denominator


def h(theta):
    return np.sin(theta) * (theta - np.sin(theta)) / (1 - np.cos(theta)) ** 2


def q(v, h, q_2):
    K = q_2 * v * h ** 2 / (1 - v + h * v) ** 2 - 2 * v * (1 - v) * (1 - h) ** 2 / (1 - v + h * v) ** 2
    return K


def q_2(theta):
    return (1 - np.cos(theta)) / np.sin(theta) ** 2


f_2 = []
T = []
Q =[[],[],[],[],[],[],[]]
V2 =[[],[],[],[],[],[],[]]
numV=7
Q_2 = []
Q_1 = []
V_2=[]
V_1=[]
colorcode=['b','g','r','c','m','y','b','g','r','c','m','y','b','g','r','c','m','y']
plt.figure()
for i in range(numV):
    f_2.append(round(0.1*(1+i),1))
    for theta in np.arange(0.00001, 5, 0.001):
        V = v(theta, f_2[i])
        H = h(theta)
        Q2 = q_2(theta)
        Q[i].append(q(V, H, Q2))
        #Q_2.append(Q2)
        #Q_1.append(0)
        V_2.append(V)
        V_1.append(1-V)
        V2[i].append(V)
        if i==0:
            T.append(theta)
i=0     
Qq=np.transpose(Q) 
for i in range(numV):
    plt.plot(T, V2[i], label=r'$v(\theta)$ for $v_t=$'+f'{f_2[i]}')
    #plt.plot(T, Q[i], label=r'$q(\theta)$ for $v_t=$'+f'{f_2[i]}')  
plt.legend()  
#plt.plot(T, Q_2, 'b',label=r'$q_2(\theta)$')
#plt.plot(T, Q_1, 'g',label=r'$q_1(\theta)$')
#plt.plot(T, V_1, 'g',label='Volume fraction in region 1')
#plt.plot(T, V_2, 'b',label='Volume fraction in region 2')
#plt.plot(T, Q, label=r'$q(\theta)$, $v_t=$'+f'{f_2}')
#plt.plot(T, Q)
plt.xlim(0, 5)
#plt.axhline(y=0, color='k', linestyle='--')
plt.xlabel(r'$\theta$')
plt.ylabel(r'v($\theta$)')
#plt.ylabel(r'q($\theta$)')

plt.title(r'Volume fraction $v$ as a'+'\n'+r'function of development angle $\theta$')
plt.show()
