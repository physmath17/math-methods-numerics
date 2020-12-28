from runge_kutta_four import *
import numpy as np
from matplotlib import pyplot as plt
from scipy.special import airy, gamma
from datetime import datetime

startTime = datetime.now()

# parameters
a = 0
bl = -15
bg = 5
Nl = 1500
Ng = 5000
xl = np.linspace(a, bl, Nl)
xg = np.linspace(a, bg, Ng)
hl = (bl - a)/Nl
hg = (bg - a)/Ng

y0 = 1/(pow(3, 2/3)*gamma(2/3))
t0 = -1/(pow(3, 1/3)*gamma(1/3))
Y = np.array([y0, t0])
solnl = np.array(Y[0])               # solution, x < 0
dsolnl = np.array(Y[1])              # derivative of solution, x < 0
solng = np.array(Y[0])               # solution, x > 0
dsolng = np.array(Y[1])              # derivative of solution, x > 0

def F(t, Y) :
    f = np.array([Y[1], t*Y[0]])
    return f

# solution in the range x < 0
for i in range(1, len(xl)) :
    Y = RK4(xl[i-1], hl, Y, F)
    solnl = np.append(solnl, Y[0])
    dsolnl = np.append(dsolnl, Y[1])

# solution in the range x > 0
Y = np.array([y0, t0])

for i in range(1, len(xg)) :
    Y = RK4(xg[i-1], hg, Y, F)
    solng = np.append(solng, Y[0])
    dsolng = np.append(dsolng, Y[1])

x = np.append(xl[::-1], xg[1::])
soln = np.append(solnl[::-1], solng[1::])

ai, aip, bi, bip = airy(x)
plt.plot(x, ai, 'r', label='Ai(x)')
plt.plot(x, soln, 'b', label='y(x)')
plt.xlabel("x")
plt.title("Airy equation solution using fourth order Runge-Kutta method")
plt.grid()
plt.legend()
plt.show()