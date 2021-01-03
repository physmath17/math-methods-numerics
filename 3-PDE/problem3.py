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

# parameters
a = 0
b = 1
hx = 0.05
ht = 0.0025

# initial conditions specfied by
phi = lambda y : ((0 < y <= 0.5 and 2*y) or 2*(1 - y))

# solving the equation
def solver(xi, xf, dx, dt, s=1, D=1) :
    ''' D : diffusion constant, xi, xf : spatial boundary points, dx : spatial increment, dt : time increment, s = dt/(dx)^2
    returns a 2D array as solution where the columns for a fixed rwo specify solutions at different times and the rows for a fixed column specify the solution at different spatial points '''
    if s != 1 :
        dt = s*dx**2/D

    Nx = int((xf-xi)/dx)
    Nt = 1000

    x = np.linspace(xi, xf, Nx)
    t = np.array([i*dt for i in range(Nt)])
    u = np.empty((Nx, Nt))

    # initial condition
    # for t theer is one initial condition is specified by phi
    for i in range(1, len(x)-1) :
        u[0][i] = phi(x[i])
    
    # for x we need to specify two initial conditions
    u[:, :1] = 0
    
    # difference equation - recurrence relation for the heat equation
    for time in t :
        for space in range(1, Nx, dx) :
            u[time][space] = s*(u[time-1][space+1] + u[time-1][space-1]) + (1 - 2*s)*u[time-1][space]
    return u

# def heatmap(v, dt, k) :
#     plt.clf()
#     plt.title(f"Temperature at t = {k*dt:.3f} unit time")
#     plt.xlabel("x")
#     plt.pcolormesh(v, cmap=plt.cm.jet, vmin=0, vmax=100)
#     plt.colorbar()

#     return plt

T = solver(a, b, hx, ht)

# def animate(k) :
#     heatmap(T[k], ht, k)

# anim = animation.FuncAnimation(plt.figure(), animate, interval=1, frames=1000, repeat=False)
# anim.save("heat_equation_solution.gif")
