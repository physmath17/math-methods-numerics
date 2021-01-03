from scipy.special import airy, gamma
import matplotlib.pyplot as plt
import numpy as np
from math import factorial as fact
from datetime import datetime

startTime = datetime.now()

v = np.linspace(-15, 15, 31)
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
            s += pow(3, -2/3)*pow(x, n)/(pow(3, 2*n/3)*fact(n/3)*gamma(n/3 + 2/3))
        elif n%3 == 1 :
            s += -pow(3, -4/3)*pow(x, n)/(pow(3, 2*(n-1)/3)*fact((n-1)/3)*gamma(n/3 + 1))
        # else:
        #     s += 0
        # print("% .4f\t\t% .4f\t%d" % (y, s, n))
        if abs(s - y) < 1e-4 :
            break
        n += 1
    return n, s

print("x\t\tn")
# approx = np.array([airy_first(v[i], ai[i]) for i in range(len(v))])
for i in range(len(v)) :
    try :
        print(v[i], "\t", airy_first(v[i], ai[i]), "\t", ai[i])
    except :
        print("Something went wrong with {}".format(v[i]))
endTime = datetime.now()
print("\nExecution time : ", endTime - startTime)

plt.scatter(v, ai)
# plt.scatter(v, approx, "r")
plt.show()