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


f_2 = 0.3
H_1 = 1
T = []
Q = []
Q_2 = []
Q_1 = []
V_2=[]
V_1=[]
for theta in np.arange(0.00001, 5, 0.001):
    V = v(theta, f_2)
    H = h(theta)
    Q2 = q_2(theta)
    Q.append(q(V, H, Q2))
    Q_2.append(Q2)
    Q_1.append(0)
    V_2.append(V)
    V_1.append(1-V)
    #Q_Q.append(H_1 * (1 - V + V * H))
    T.append(theta)
#plt.plot(T, Q, 'b',label=r'$q(\theta)$')
#plt.plot(T, Q_2, 'b',label=r'$q_2(\theta)$')
#plt.plot(T, Q_1, 'g',label=r'$q_1(\theta)$')
plt.plot(T, V_1, 'g',label='Volume fraction in region 1')
plt.plot(T, V_2, 'b',label='Volume fraction in region 2')

plt.xlim(0, 5)
#plt.axhline(y=0, color='r', linestyle='--')
plt.xlabel(r'$\theta$')
plt.ylabel(r'v($\theta$)')
#plt.ylabel(r'q($\theta$)')
plt.legend()
plt.title(r'Volume fraction v as a'+'\n'+r'function of development angle $\theta$')
plt.show()
