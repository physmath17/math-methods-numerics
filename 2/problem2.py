'''
First part of the problem - solve the Airy equation and plot the results
Second part of the problem - use series soulution to approximate the Airy function at x = -10000, -100, -1, 1, 100, 10000
'''

from scipy.special import airy, gamma
import matplotlib.pyplot as plt
import numpy as np
from math import factorial as fact
from datetime import datetime

startTime = datetime.now()

# parameters
a = 0 
# y0 = 1/(pow(3, 2/3)*gamma(2/3))         # for Ai(x)
# dy0 = -1/(pow(3, 1/3)*gamma(1/3))       # for Ai(x)
y0 = 1/(pow(3, 1/6)*gamma(2/3))         # for Bi(x)
dy0 = pow(3, 1/6)/gamma(1/3)            # for Bi(x)

# solving the Airy equation using central difference
"""
def cntrl(xi, xf, yi, dyi, steps) :
    ''' solves the Airy equation y'' = xy
    xi : initial x datapoint, xf : final x datapoint, yi : initial y datapoint, dyi : derivativative initial data, steps : number of iterations,
    returns two 1D array of x, y values '''

    h = (xf - xi)/steps
    x = np.linspace(xi, xf, steps)
    y = np.array([yi])

    # using disrete derivative implementation of the continuous derivative
    y1 = (h*dyi + yi)
    y = np.append(y, y1)
    for i in range(0, len(x)-2) :
        dy = (2 + x[i + 1]*h**2)*y[i + 1] - y[i]
        y = np.append(y, dy)

    return x, y

# solving the Airy equation using forward difference
def frwrd(xi, xf, yi, dyi, steps) :
    ''' solves the Airy equation y'' = xy
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

xg, yg = cntrl(a, 2, y0, dy0, 1000)      # x range greater than 0
xl, yl = cntrl(a, -15, y0, dy0, 1500)    # x range less than 0

x = np.append(xl[::-1], xg[1::])
y = np.append(yl[::-1], yg[1::])
"""
v = [-10000, -100, -1, 1, 100, 10000]
ai, aip, bi, bip = airy(v)

# approximating the Airy equation with the series solution
def airy_first(x, y) :
    ''' returns the number of terms required in the power series to approximate the Airy function of the first kind at x and and the value of Ai(x) '''

    s = 0.0
    n = 0
    while True :
        if n%3 == 0 :
            s += pow(3, -2/3)*pow(x, n)/(pow(3, 2*n/3)*fact(n/3)*gamma(n/3 + 2/3))  
        elif n%3 == 1 :
            s += -pow(3, -4/3)*pow(x, n)/(pow(3, 2*(n-1)/3)*fact((n-1)/3)*gamma(n/3 + 1))
        else:
            s += 0

        if abs(s - y) < 1e-3 :
            break
        n += 1
    return n, s

def airy_second(x, y) :
    ''' returns the number of terms required in the power series to approximate the Airy function of the first kind at x and and the value of Ai(x) '''

    s = 0.0
    n = 0
    while True :
        s += pow(3, -1/6)*pow(x, n)/(pow(3, 2*n/3)*fact(n/3)*gamma(n/3 + 2/3))*(n%3 == 0) - pow(3, -5/6)*pow(x, n)/(pow(3, 2*(n-1)/3)*fact((n-1)/3)*gamma(n/3 + 1))*(n%3 == 1) + 0*(n%3 == 2)

        if abs(s - y) < 1e-3 :
            break
        n += 1
    return n, s
x=1
s=0
for n in range(600) :
    if n%3 == 0 :
            s += pow(3, -2/3)*pow(x, n)/(pow(3, 2*n/3)*fact(n/3)*gamma(n/3 + 2/3))  
    elif n%3 == 1 :
        s += -pow(3, -4/3)*pow(x, n)/(pow(3, 2*(n-1)/3)*fact((n-1)/3)*gamma(n/3 + 1))
    else:
        s += 0
print(s, ai[4])

endTime = datetime.now()
print("Execution time : ", endTime - startTime)

# plt.plot(x, ai, 'r', label='Ai(x)')
# plt.plot(x, bi, 'r', label='Bi(x)')
# plt.plot(x, y, 'b', label='y(x)')
# plt.ylim(-0.5, 0.6)             # for Ai(x)
# plt.ylim(-0.6, 4)                 # for Bi(x)
# plt.xlabel("x")
# plt.title("Airy equation solution using central difference method")
# plt.grid()
# plt.legend(loc='upper left')
# plt.show()