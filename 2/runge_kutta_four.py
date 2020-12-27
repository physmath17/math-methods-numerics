"""
THe Runge-Kutta method for a second order ODE is solved by 
converting the problem into two linear ODEs. 

In general if we have m linear ODEs:

y1' = f1(x,y1,y2,...,ym)
y2' = f2(x,y1,y2,...,ym)
.
.
.
ym' = fm(x,y1,y2,...,ym)

define Yvec = (y1, y2, ..., ym) and Fvec = (f1, f2, ..., fm)

So the systen of equations can be written as:

Yvec' = Fvec(x, Yvec)

The genraliezed RK4 is then given by:

k1vec = h*Fvec(x_n, Yvec(x_n))
k2vec = h*Fvec(x_n + h/2, Yvec(x_n) + k1vec/2)
k3vec = h*Fvec(x_n + h/2, Yvec(x_n) + k2vec/2)
k4vec = h*Fvec(x_n + h, Yvec(x_n) + k3vec)

The solution is then given by:

Yvec(x_n+1) = Yvec(x_n) + (k1vec + 2*k2vec + 2*k3vec + k4vec)/6

with m initial conditiond specified by Yvec(x_0)

"""

import numpy as np

def RK4(x, h, Yvec, Fvec):
   k1 = h*Fvec(x, Yvec)
   k2 = h*Fvec(x + h/2, Yvec + k1/2)
   k3 = h*Fvec(x + h/2, Yvec + k2/2)
   k4 = h*Fvec(x + h, Yvec + k3)
    
   return (Yvec + (k1 + 2*k2 + 2*k3 + k4)/6)