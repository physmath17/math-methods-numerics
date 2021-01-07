'''
solving the advection equation in one dimension for Dirichlet boundary conditions ---- (take D = 1)
u_t + v*u_x = 0, with -L < x < L, t > 0,  
u(x, 0) = phi(x) = 0       if -L < x < 0
                 = 0.5     if x = 0
                 = 1       if 0 < x < L
'''

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.animation import FuncAnimation
from mpl_toolkits.mplot3d import Axes3D

# parameters
L = 10
hx = 0.05
ht = 0.05
v = 1

# initial conditions specfied by
phi = lambda y : 0*(y<0) + 0.5*(y==0) + 1*(y>0)

# solving the equation
def advection_solver(phi, l, dx, dt, v=1) :
    ''' phi : specifies the boundary conditions, l : boundary point, dx : spatial increment, dt : time increment, v = velocity
    returns a 2D array as the solution to the advection equation where the rows for a fixed column specify solutions at different times and the columns for a fixed row specify the solution at different spatial points '''
    
    s = v*dt/(2*dx)

    Nx = int(2*l/dx) + 1
    Nt = 201

    x = np.linspace(-l, l, Nx)
    t = np.array([i*dt for i in range(Nt)])
    u = np.empty((Nt, Nx))

    # initial condition
    # for t there is one initial condition is specified by phi
    for i in range(Nx) :
        u[0][i] = phi(x[i])
    
    for i in range(Nt) :
        u[i][0] = 0
        u[i][Nx-1] = 1

    # difference equation - recurrence relation for the heat equation
    for time in range(Nt-1) :
        for space in range(1, Nx-1) :
            u[time+1][space] = u[time][space] - s*(u[time][space+1] + u[time][space-1])
    return x, t, u

def modified_advection_solver(phi, l, dx, dt, v=1) :
    ''' phi : specifies the boundary conditions, l : boundary point, dx : spatial increment, dt : time increment, v = velocity
    returns a 2D array as the solution to the advection equation where the rows for a fixed column specify solutions at different times and the columns for a fixed row specify the solution at different spatial points '''
    
    D = v*dx/2
    r = v*dt/(2*dx)
    s = dt*D/(dx**2)

    Nx = int(2*l/dx) + 1
    Nt = 201

    x = np.linspace(-l, l, Nx)
    t = np.array([i*dt for i in range(Nt)])
    u = np.empty((Nt, Nx))

    # initial condition
    # for t there is one initial condition is specified by phi
    for i in range(Nx) :
        u[0][i] = phi(x[i])
    
    for i in range(Nt) :
        u[i][0] = 0
        u[i][Nx-1] = 1

    # difference equation - recurrence relation for the heat equation
    for time in range(Nt-1) :
        for space in range(1, Nx-1) :
            u[time+1][space] = u[time][space] - r*(u[time][space+1] - u[time][space-1]) + s*(u[time][space+1] - 2*u[time][space] + u[time][space-1]) 
    return x, t, u

x, t, u = modified_advection_solver(phi, L, hx, ht, v)

def field(x, f, k) :
    plt.clf()
    plt.title(r"Field at t = {time:.3f} unit time for $\Delta$t = {dt}, $\Delta$x = {dx},  v = {vel}".format(time=ht*k, dt=ht, dx=hx, vel=v)+"\n")
    plt.xlabel("x")
    plt.ylabel("u")
    plt.plot(x, f)
    return plt

def animate(k) :
    return field(x, u[k], k)

fig, ax = plt.subplots()
anim = animation.FuncAnimation(fig, animate, interval=200, frames=200, repeat=False)
anim.save("advection_equation_solution_step-function_v_1.gif", writer='imagemagick', fps=10)

# fig = plt.figure() 
# ax = fig.add_subplot(111, projection='3d')
# space, time = np.meshgrid(x, t)
# # ax.plot3D(space, time, u[0], '-k')
# ax.plot_surface(space, time, u)
# ax.set_title(r"Field at different positions and times for $\Delta$t = {dt}, $\Delta$x = {dx},  v = {vel}".format(dt=ht, dx=hx, vel=v)+"\n")
# ax.set_xlabel('x')
# ax.set_ylabel('t')
# ax.set_zlabel('field amplitude')
# plt.show()