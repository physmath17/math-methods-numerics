import numpy as np
import matplotlib.pyplot as plt

for i in range(10, 61, 10) :
    data = np.loadtxt("data_L_{}.txt".format(i))
    plt.scatter(data[:,0], i**0.1*abs(data[:,2]))
plt.show()