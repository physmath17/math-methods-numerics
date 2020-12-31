from scipy.special import airy, gamma
import matplotlib.pyplot as plt
import numpy as np
from math import factorial as fact
from datetime import datetime

startTime = datetime.now()

v = [-10000, -100, -10, -1, 0, 1, 10, 100, 10000]
ai, aip, bi, bip = airy(v)

s = 0
x = 1
for n in range(10) :
    if n%3 == 0 :
            s += pow(3, -2/3)*pow(x, n)/(pow(9, n/3)*fact(n/3)*gamma(n/3 + 2/3))  
    elif n%3 == 1 :
        s += -pow(3, -4/3)*pow(x, n)/(pow(9, (n-1)/3)*fact((n-1)/3)*gamma(n/3 + 1))
    elif n%3 == 2 :
        s += 0

    print("%f\t%f" % (s, ai[5]))

endTime = datetime.now()
print("\nExecution time : ", endTime - startTime)
