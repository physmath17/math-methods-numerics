'''
First part of the problem - solve the Airy equation and plot the results
Second part of the problem - use series soulution to approximate the Airy function and plot the approximation along with the original function
'''

from scipy.special import airy, gamma
import matplotlib.pyplot as plt
import numpy as np
from math import exp

# parameters
a = 0 
y0 = 1/(pow(3, 2/3)*gamma(2/3))
dy0 = -1/(pow(3, 1/3)*gamma(1/3))

# solving the Airy equation using central difference
def cntrl(xi, xf, yi, dyi, steps) :
    ''' solves the Airy equation y''(x) = xy
    xi : initial x datapoint, xf : final x datapoint, yi : initial y datapoint, dyi : derivativative initial data, steps : number of iterations,
    returns two 1D array of x, y values '''

    h = (xf - xi)/steps
    x = np.linspace(xi, xf, steps)
    y = np.array([yi])

    # using disrete derivative implementation of the continuous derivative
    y1 = (h*dyi + yi)
    y = np.append(y, y1)
    for i in range(0, len(x)-2) :
        dy = (2 + x[i]*h**2)*y[i + 1] - y[i]
        y = np.append(y, dy)

    return x, y

# solving the Airy equation using forward difference
def frwrd(xi, xf, yi, dyi, steps) :
    ''' solves the Airy equation y''(x) = xy
    xi : initial x datapoint, xf : final x datapoint, yi : initial y datapoint, dyi : derivativative initial data, steps : number of iterations,
    returns two 1D array of x, y values '''

    h = (xf - xi)/steps
    x = np.linspace(xi, xf, steps)
    y = np.array([yi])

    # using disrete derivative implementation of the continuous derivative
    y1 = (h*dyi + yi)
    y = np.append(y, y1)
    for i in range(0, len(x)-2) :
        dy = (x[i]*h**2 - 1)*y[i] + 2*y[i + 1]
        y = np.append(y, dy)

    return x, y

xg, yg = cntrl(a, 10, y0, dy0, 1000)       # x range greater than 0
xl, yl = cntrl(a, -15, y0, dy0, 1500)    # x range less than 0

x = np.append(xl[::-1], xg[1::])
y = np.append(yl[::-1], yg[1::])
ai, aip, bi, bip = airy(x)

# approximating the Airy equation with the series solution

print(y[-1])

plt.plot(x, ai, 'k', label='Ai(x)')
plt.plot(x, y, 'b', label='y(x)', alpha=0.5)
# plt.plot(x, bi, 'b--', label='Bi(x)')
plt.ylim(-0.5, 1.0)
plt.grid()
plt.legend(loc='upper left')
plt.show()