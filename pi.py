import math
from mpmath import mp

mp.dps = 20

def factorial(n) :
    """ compute n! """
    if n == 0 :
        return 1
    else :
        return n*factorial(n-1)

def estimate_pi() :
    """ uses the infinte series given by Ramanujan to estimate the value of pi """
    s = 0.0
    k = 0
    while True :
        s = s + factorial(4*k)*(1103 + 26390*k) / (pow(factorial(k),4) * pow(396, 4*k))
        pi = 9801 / (2*s*math.sqrt(2))
        
        if abs(pi - math.pi) < 1e-21 :
            break
        k += 1
    print("The value of pi is (using Ramanujan's series) = {}".format(mp.pi))
    print("The value of pi is (using Gaussian integral) = {}".format(mp.quad(lambda x: mp.exp(-x**2), [-mp.inf, mp.inf]) ** 2))

estimate_pi()

