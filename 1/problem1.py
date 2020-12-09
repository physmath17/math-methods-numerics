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
N = int(input("Enter the number of steps : "))
a = 0 #float(input("Enter the intial point : "))
b = 2*np.pi #float(input("Enter the final point : "))
y0 = 0 #float(input("Enter the value of function at the initial point : "))
dy0 = 1 #float(input("Enter the value of the derivative of the function at the initial point : "))
h = (b - a)/N
x = np.linspace(a, b, N)
y = np.array([y0])

# actual solution
sin = np.array([np.sin(x[i]) for i in range(len(x))])

# usinf forward derivative implementation of the continuous derivative
y1 = (h*dy0 + y0)
y = np.append(y, y1)

for i in range(0, len(x)-2) :
    dy = 2*y[i + 1] - (1 + h**2)*y[i]
    y = np.append(y, dy)

plt.plot(x, sin, color='black')
plt.scatter(x, y, color='red', s=2)
plt.title("Solution to y'' = -y, y(0)=0, y'(0)=1 using forward discrete derivative with N = 1000")
plt.xlabel("x")
plt.legend(["sin x", "y(x)"])
plt.show()