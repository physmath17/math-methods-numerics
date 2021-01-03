'''
solving the heat equation in one dimension for Dirichlet boundary conditions ---- (take D = 1)
u_t - D*u_xx = 0, with 0 < x < 1, t > 0, 
u(0, t) = 0 = u(1, t), 
u(x, 0) = phi(x) = 2*x       if 0 < x <= 0.5
                 = 2*(1 - x) if 0.5 <= x < 1.
'''

import numpy as np
import matplotlib.pyplot as plt

# parameters
a = 0
b = 1
hx = 0.05
ht = 0.0025

# initial conditions specfied by
phi = lambda y : ((0 < y <= 0.5 and 2*y) or 2*(1 - y))

# solving the equation
def solver(D=1, xi, xf, dx=0.05, dt=0.0025, s=1) :
    ''' D : diffusion constant, xi, xf : spatial boundary points, dx : spatial increment, dt : time increment, s = dt/(dx)^2
    returns a 2D array as solution where the columns for a fixed rwo specify solutions at different times and the rows for a fixed column specify the solution at different spatial points '''
    if s != 1 :
        dt = s*dx**2

    Nx = int((b-a)/dx)
    Nt = 1000

    x = np.linspace(0, 1, Nx)
    t = np.array([i*dt for i in range(Nt)])
    u = np.empty((Nx, Nt))

    # initial condition
    # for t theer is one initial condition is specified by phi
    for i in range(1, len(x)-1) :
        u[i][0] = phi(x[i])
    
    # for x we need to specify two initial conditions
    for i in range(len(t)) :
        u[0][i] = 0
        u[-1][i] = 0 
    
    # difference equation - recurrence relation for the heat equation
    for time in range(1, Nt) :
        for space in range(1, Nx) :
            u[space][time] = s*(u[space+1][time-1] + u[space-1][time-1]) + (1 - 2*s)*u[space][time-1]
    
    fig, ax = plt.subplots()