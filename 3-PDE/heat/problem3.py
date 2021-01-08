'''
solving the heat equation in one dimension for Dirichlet boundary conditions ---- (take D = 1)
u_t - D*u_xx = 0, with 0 < x < 1, t > 0, 
u(0, t) = 0 = u(1, t), 
u(x, 0) = phi(x) = 2*x       if 0 < x <= 0.5
                 = 2*(1 - x) if 0.5 <= x < 1.
'''

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.animation import FuncAnimation
from mpl_toolkits.mplot3d import Axes3D

# parameters
a = 0
b = 1
hx = 0.05
ht = 0.0025

# initial conditions specfied by
phi = lambda y : 2*y*(y<0.5) + 2*(1-y)*(y>=0.5)

# solving the equation
def heat_solver(phi, xi, xf, dx=0.05, dt=0.0025, s=1, D=1) :
    ''' phi : specifies the boundary conditions, xi, xf : spatial boundary points, dx : spatial increment, dt : time increment, s = D*dt/(dx)^2, D : diffusion constant
    returns a 2D array as the solution to the heat equation where the rows for a fixed column specify solutions at different times and the columns for a fixed row specify the solution at different spatial points '''
    if s != 1 and dx == 0.05 :
        dt = s*(dx**2)/D
    elif s!= 1 and dt == 0.0025 :
        dx = np.sqrt(D*dt/s)
    else :
        s = D*dt/(dx**2)

    Nx = int((xf-xi)/dx) + 1
    Nt = 301

    x = np.linspace(xi, xf, Nx)
    t = np.array([i*dt for i in range(Nt)])
    u = np.empty((Nt, Nx))

    # initial condition
    # for t there is one initial condition is specified by phi
    for i in range(Nx) :
        u[0][i] = phi(x[i])
    
    # for x we need to specify two initial conditions
    for i in range(1, Nt) :
        u[i][0] = 0
        u[i][Nx-1] = 0

    # difference equation - recurrence relation for the heat equation
    for time in range(1, Nt) :
        for space in range(1, Nx-1) :
            u[time][space] = s*(u[time-1][space+1] + u[time-1][space-1]) + (1 - 2*s)*u[time-1][space]
    return x, t, u, dt

def heatmap(x, v, dt, k) :
    plt.clf()
    plt.title(f"Temperature at t = {k*dt:.3f} unit time")
    plt.xlabel("x")
    plt.ylabel("u")
    plt.plot(x, v)
    return plt

x, t, u, delta_t = heat_solver(phi, a, b, dt=ht, s=0.48)

def animate(k) :
    return heatmap(x, u[k], delta_t, k)

# fig, ax = plt.subplots()
# anim = animation.FuncAnimation(fig, animate, interval=200, frames=300, repeat=False)
# anim.save("heat_equation_solution_s_0.51.gif", writer='imagemagick', fps=10)

fig = plt.figure() 
ax = fig.add_subplot(111, projection='3d')
space, time = np.meshgrid(x, t) 
ax.plot_surface(space, time, u)
ax.set_title((f"Temperature distribution in a rod\n"))
ax.set_xlabel('x')
ax.set_ylabel('t')
ax.set_zlabel('temperature')
plt.show()