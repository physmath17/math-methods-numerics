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
"""
# parameters
a = 0 
y0 = 1/(pow(3, 2/3)*gamma(2/3))         # for Ai(x)
dy0 = -1/(pow(3, 1/3)*gamma(1/3))       # for Ai(x)
# y0 = 1/(pow(3, 1/6)*gamma(2/3))         # for Bi(x)
# dy0 = pow(3, 1/6)/gamma(1/3)            # for Bi(x)

# solving the Airy equation using central difference
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

xg, yg = cntrl(a, 6, y0, dy0, 100000)      # x range greater than 0
xl, yl = cntrl(a, -15, y0, dy0, 1500)    # x range less than 0

x = np.append(xl[::-1], xg[1::])
y = np.append(yl[::-1], yg[1::])
ai, aip, bi, bip = airy(x)
"""
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

        if abs(s - y) < 1e-4 :
            break
        n += 1
    return s, n

def airy_second(x, y) :
    ''' returns the number of terms required in the print("%f\t%f" % airy_first(1, ai[3]), ai[3])
power series to approximate the Airy function of the first kind at x and and the value of Ai(x) '''

    s = 0.0
    n = 0
    while True :
        if n%3 == 0 :
            s += pow(3, -1/6)*pow(x, n)/(pow(3, 2*n/3)*fact(n/3)*gamma(n/3 + 2/3))  
        elif n%3 == 1 :
            s += -pow(3, -5/6)*pow(x, n)/(pow(3, 2*(n-1)/3)*fact((n-1)/3)*gamma(n/3 + 1))
        else:
            s += 0

        if abs(s - y) < 1e-4 :
            break
        n += 1
    return s, n

v = [-10000, -100, -10, -1, 0, 1, 10, 100, 10000]
ai, aip, bi, bip = airy(v)

print("x\t Ai(x)\t\t approximation\t steps")
for i in range(2, 7) :
    approx, terms_reqd = airy_first(v[i], ai[i])
#    print("{}\t {:.5f}\t  {:.5f}\t {}".format(v[i], ai[i], approx, terms_reqd))
    print("% d\t% .5f\t% .5f\t %d" % (v[i], ai[i], approx, terms_reqd))

# print("\nx\t Bi(x)\t\t approximation\t steps")
# for i in range(3, len(v)-3) :
#     approx, terms_reqd = airy_second(v[i], bi[i])
#     print("% d\t% .5f\t% .5f\t %d" % (v[i], bi[i], approx, terms_reqd))
print(airy_second(1, bi[5]))

endTime = datetime.now()
print("\nExecution time : ", endTime - startTime)

# plt.plot(x, ai, 'r', label='Ai(x)')
# plt.plot(x, bi, 'r', label='Bi(x)')
# plt.plot(x, y, 'b', label='y(x)')
# plt.ylim(-0.5, 0.6)             # for Ai(x)
# plt.ylim(-0.6, 4)               # for Bi(x)
# plt.xlabel("x")
# plt.title("Airy equation solution using central difference method")
# plt.grid()
# plt.legend(loc='upper left')
# plt.show()