"""
Solving y'' = -y, y(0) = 0, y'(0) = 1

General scheme :

y'' = g ---------> dx_i = x_i - x_(i-1) = dx = h, g_n = (y_(n+1) - 2y_n + y_(n-1))/h^2 --> y_n = (y_(n+1) + y_(n-1) - g_n*h^2)/2


forward derivative : dfy = (y(x + dx) - y(x))/dx ----> (y_(n+1) - y_n)/dx
backward derivative : dby = (y(x) - y(x - dx))/dx ----> (y_(n) - y_(n-1))/dx
discrete derivatve : ddy = (y(x + dx) - y(x - dx))/2dx
"""

import numpy as np
import matplotlib.pyplot as plt

# parameters
# N = int(input("Enter the number of steps : "))
a = 0 #float(input("Enter the intial point : "))
b = 2*np.pi #float(input("Enter the final point : "))
y0 = 0 #float(input("Enter the value of function at the initial point : "))
dy0 = 1 #float(input("Enter the value of the derivative of the function at the initial point : "))

def frwrd(xi, xf, yi, dyi, steps) :
    ''' xi : initial x datapoint, xf : final x datapoint, yi : initial y datapoint, dyi : derivativative initial data, steps : number of iterations,
    returns two 1D array of x, y values '''

    h = (xf - xi)/steps
    x = np.linspace(xi, xf, steps)
    y = np.array([yi])

    # using forward difference implementation of the continuous derivative
    y1 = (h*dyi + yi)
    y = np.append(y, y1)
    for i in range(0, len(x)-2) :
        dy = 2*y[i + 1] - (1 + h**2)*y[i]
        y = np.append(y, dy)

    return x, y

def dscrt(xi, xf, yi, dyi, steps) :
    ''' xi : initial x datapoint, xf : final x datapoint, yi : initial y datapoint, dyi : derivativative initial data, steps : number of iterations,
    returns two 1D array of x, y values '''

    h = (xf - xi)/steps
    x = np.linspace(xi, xf, steps)
    y = np.array([yi])

    # using central difference implementation of the continuous derivative
    y1 = (h*dyi + yi)
    y = np.append(y, y1)
    for i in range(0, len(x)-2) :
        dy = (2 - h**2)*y[i + 1] - y[i]
        y = np.append(y ,dy)

    return x, y

# results
step = [1000, 2000, 5000, 10000]
for N in step :
    x, y = frwrd(a, b, y0, dy0, N)
    plt.scatter(x, y, s=2)

# actual solution
    u = np.linspace(a, b, 1000)
    sin = np.array([np.sin(u[i]) for i in range(len(u))])

plt.plot(u, sin, color='black')
plt.title("Solution to y'' = -y, y(0)=0, y'(0)=1 using forward difference")
# plt.title("Solution to y'' = -y, y(0)=0, y'(0)=1 using central difference")
plt.xlabel("x")
plt.legend(["sin x", "y(x) with N = 1000", "y(x) with N = 2000", "y(x) with N = 5000", "y(x) with N = 10000"], loc='upper right')
plt.show()