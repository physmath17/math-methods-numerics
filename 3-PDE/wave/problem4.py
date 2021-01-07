'''
solving the wave equation in one dimension for Dirichlet boundary condition in space and mixed boundary condition in time ---- (take c = 1)
u_tt - c^2*u_xx = 0, with 0 < x < L, t > 0, 
u(0, t) = 0 = u(L, t), 
u(x, 0) = phi(x) 
u_t(x, 0) = psi(x)
'''

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.animation import FuncAnimation
from mpl_toolkits.mplot3d import Axes3D

# parameters
a = 0
L = 20
Nx = 2001
Nt = 200

# initial conditions specfied by
phi = lambda y : np.exp(-(y - L/2)**2) - np.exp(-(L/2)**2)
psi = lambda y : 0

# solving the equation
def wave_solver(phi, psi, xi, xf, Nx, Nt, s) :
    ''' phi : specfies u(x, 0), psi : specifies u_t(x, 0), xi, xf : spatial boundary points, N : number of spatial intervals, s : square of ratio of time interval to spatial interval (can be used to define the time interval), c : speed of wave, by default set to ubnity
    returns a 2D array as the solution to the wave equation where the rows for a fixed column specify solutions at different times and the columns for a fixed row specify the solution at different spatial points '''

    dx = (xf - xi)/Nx
    dt = dx*s**0.5

    x = np.linspace(xi, xf, Nx)
    t = np.array([i*dt for i in range(Nt)])
    u = np.empty((Nt, Nx))

    # initial condition
    # for t there are two initial conditions specified by phi and psi
    for i in range(1, Nx-1) :
        u[0][i] = phi(x[i])
        u[1][i] = psi(x[i])*dt + (1 - s)*phi(x[i]) + 0.5*s*(phi(x[i+1]) + phi(x[i-1]))
    
    # for x we need to specify two initial conditions
    for i in range(Nt) :
        u[i][0] = 0
        u[i][Nx-1] = 0

    # difference equation - recurrence relation for the wave equation
    for time in range(1, Nt-2) :
        for space in range(1, Nx-1) :
            u[time+1][space] = s*(u[time][space+1] + u[time][space-1]) + 2*(1 - s)*u[time][space] - u[time-1][space]
    return x, t, u, dt

def wave(x, v, dt, k) :
    plt.clf()
    plt.title(f"Wave amplitude at t = {k*dt:.3f} unit time")
    plt.xlabel("x")
    plt.ylabel("u")
    plt.plot(x, v)
    return plt

x, t, u, delta_t = wave_solver(phi, psi, a, L, Nx, Nt, 0.9)
def animate(k) :
    return wave(x, u[k], delta_t, k)

# fig, ax = plt.subplots()
# anim = animation.FuncAnimation(fig, animate, interval=200, frames=200, repeat=False)
# anim.save("wave_equation_solution_s_1.1.gif", writer='imagemagick', fps=10)

fig = plt.figure() 
ax = fig.add_subplot(111, projection='3d')
space, time = np.meshgrid(x, t) 
ax.plot_surface(space, time, u)
ax.set_title((f"Wave solution at different positions and times\n"))
ax.set_xlabel('x')
ax.set_ylabel('t')
ax.set_zlabel('amplitude')
plt.show()