import scipy as sc
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import numpy as np


def r_move(rx, vx, dt):
    nrx = rx + vx * dt
    return nrx


def v_move(vx, ax, dt):
    nvx = vx + ax * dt
    return nvx


height = 30  # meters
length = 13.7  # meters 45 feet
v_max = 10  # m/s
dt = 0.001
t_tot = 10
length = int(t_tot / dt)
theta = np.pi / 4

# the initial condit
r = np.zeros((2, length))
v = np.zeros((2, length))
r[1][0] = height
v[0][0] = v_max * np.cos(theta)
v[1][0] = v_max * np.sin(theta)
a_x = 0
a_y = -9.82  # m/s**2

a = np.zeros((2, 1))
a[1] = a_y
t_end = 0
for t in range(round(t_tot / dt)):
    for dirc in range(2):
        r[dirc][t + 1] = r_move(rx=r[dirc][t], vx=v[dirc][t], dt=dt)
        v[dirc][t + 1] = v_move(vx=v[dirc][t], ax=a[dirc], dt=dt)

    if r[1][t] <= 0:
        t_end = t
        print(f'xend ={r[0][t]} and yend = {r[1][t]}')
        break


def dSdt2(S,t,ax,ay):
    '''
    The important part is that the vector S is defined so that it is the lists we
    want in the end, like the x pos and the y pos and the velocity in the different
    directions fx.
    '''
    x,y,vx,vy = S
    '''
    then after that we define de differential equations that we have, we
    do not have the euler equations, it is just the pure diff eq's with nothing extra to them
    '''
    dxdt = vx
    dydt = vy
    dvxdt = ax
    dvydt = ay

    '''
    then we just return the values of those diff eq's in the list, this is done many times, which
    is also why the ans has the formatting it has.
    '''
    return [dxdt,dydt,dvxdt,dvydt]


time_tot = np.linspace(0, t_tot, length)
init_condt = [r[0][0], r[1][0], v[0][0], v[1][0]]
args_list = (a[0], a[1])

ans = odeint(dSdt2, y0=init_condt, t=time_tot, args=args_list)
ans = ans.T
rx = ans[0]
ry = ans[1]

plt.figure()
plt.plot(r[0][0:t_end], r[1][0:t_end], 'r*', label='Euler integration')
plt.plot(rx, ry, 'b', label='Scipy integrator')

plt.ylim(r[1][t_end],height*1.1)
plt.xlim(0,r[0][t_end]*1.2)
plt.legend()
plt.show()
