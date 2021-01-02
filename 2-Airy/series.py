from scipy.special import airy, gamma
import matplotlib.pyplot as plt
import numpy as np
from math import factorial as fact
from datetime import datetime

startTime = datetime.now()

v = [-10000, -100, -10, -1, 0, 1, 10, 20, 100, 10000]
ai, aip, bi, bip = airy(v)

s = 0
# x = 1
# for n in range(10) :
#     if n%3 == 0 :
#             s += pow(3, -2/3)*pow(x, n)/(pow(9, n/3)*fact(n/3)*gamma(n/3 + 2/3))  
#     elif n%3 == 1 :
#         s += -pow(3, -4/3)*pow(x, n)/(pow(9, (n-1)/3)*fact((n-1)/3)*gamma(n/3 + 1))
#     elif n%3 == 2 :
#         s += 0

#     print("%f\t%f" % (s, ai[5]))

# approximating the Airy equation with the series solution
def airy_first(x, y) :
    ''' returns the number of terms required in the power series to approximate the Airy function of the first kind at x and and the value of Ai(x) '''

    s = 0.0
    n = 0
    while True :
        # t = s
        if n%3 == 0 :
            s += round(pow(3, -2/3)*pow(x, n)/(pow(3, 2*n/3)*fact(n/3)*gamma(n/3 + 2/3)),4) 
            n += 1
        elif n%3 == 1 :
            s += round(-pow(3, -4/3)*pow(x, n)/(pow(3, 2*(n-1)/3)*fact((n-1)/3)*gamma(n/3 + 1)),4)
            n += 2
        # else:
        #     s += 0
        print("% .4f\t\t% .4f\t%d" % (y, s, n))
        if abs(s - y) < 1e-4 :
            break
        # n += 1
    return s, n

# print("x\t Ai(x)\t\t approximation\t steps")
# for i in range(2, 7) :
#     approx, terms_reqd = airy_first(v[i])
#     print("% d\t% .5f\t% .5f\t %d" % (v[i], ai[i], approx, terms_reqd))
# print('\n', airy_first(15))
airy_first(20, ai[-3])
endTime = datetime.now()
print("\nExecution time : ", endTime - startTime)
