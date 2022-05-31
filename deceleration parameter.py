import numpy as np
import matplotlib.pyplot as plt
import scipy as sp

def v(theta, v_t):
    numerator = v_t * np.pi ** 3 * (1 - np.cos(theta)) ** 3
    denominator = 8 * (1 - v_t) * (theta - np.sin(theta)) ** 3 + v_t * np.pi ** 3 * (1 - np.cos(theta)) ** 3
    return numerator / denominator

def h(theta):
    return np.sin(theta) * (theta - np.sin(theta)) / (1 - np.cos(theta)) ** 2

def q(v, h, q_2):
    K = q_2 * v * h ** 2 / (1 - v + h * v) ** 2 - 2 * v * (1 - v) * (1 - h) ** 2 / (1 - v + h * v) ** 2
    return K

def q_2(theta):
    return (1 - np.cos(theta)) / np.sin(theta) ** 2

v_t = []
T = []
Q =[[],[],[],[],[],[],[]]
V2 =[[],[],[],[],[],[],[]]
numV=7
V_1=[]
plt.figure()
for i in range(numV):
    v_t.append(round(0.1*(1+i),1))
    for theta in np.arange(0.00001, 5, 0.001):
        V = v(theta, v_t[i])
        H = h(theta)
        Q2 = q_2(theta)
        Q[i].append(q(V, H, Q2))
        V_1.append(1-V)
        V2[i].append(V)
        if i==0:
            T.append(theta)
for i in range(numV):
    plt.plot(T, V2[i], label=r'$v(\theta)$ for $v_t=$'+f'{v_t[i]}')
    #plt.plot(T, Q[i], label=r'$q(\theta)$ for $v_t=$'+f'{v_t[i]}')  
#plt.plot(T, V_1, 'g',label='Volume fraction in region 1')
#plt.plot(T, Q, label=r'$q(\theta)$, $v_t=$'+f'{v_t}')
plt.legend()  
plt.xlim(0, 5)
#plt.axhline(y=0, color='k', linestyle='--')
plt.xlabel(r'$\theta$')
plt.ylabel(r'v($\theta$)')
#plt.ylabel(r'q($\theta$)')
plt.title(r'Volume fraction $v$ as a'+'\n'+r'function of development angle $\theta$')
plt.show()