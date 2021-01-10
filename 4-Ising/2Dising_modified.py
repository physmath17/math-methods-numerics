# simulating the 2D Ising model

import numpy as np
import matplotlib.pyplot as plt
from random import uniform as ran , randint, choice
from numba import jit

# parameters
L = int(input())
vol = L*L
step_size = 0.01
N = round((0.6 - 0.2)/step_size)
betaJ = np.linspace(0.2, 0.6, N+1, dtype=np.float64)  # stores inverse temperature values
E = np.zeros(N+1, dtype=np.float64)   # stores average energy
M = np.zeros(N+1, dtype=np.float64)   # stores average magnetisation
C = np.zeros(N+1, dtype=np.float64)   # stores specific heat
chi = np.zeros(N+1, dtype=np.float64) # stores susceptibility
MCsteps = 20000
n1, n2 = 1/(MCsteps*vol), 1/(MCsteps**2 * vol)

@jit(nopython=True)
def nbhr(k) :
    ''' k : index value of 1D array, dirn : 0 or 1 for horizontal and vertical direction esectively 
    returns the nearest neighbours
    '''
    i = k//L
    j = k%L
    up = (((i - 1)%L)*L + j)%vol
    down = (((i + 1)*L) + j)%vol
    right = (i*L + (j + 1)%L)%vol
    left = (i*L + (j - 1)%L)%vol
    
    return np.array([up, down, right, left])

# initial configuration
# hot start
S = np.array([choice([-1, 1]) for k in range(vol)])
# cold start
# S = np.ones(N)

@jit(nopython=True)
def Ising(S, beta, steps) : 
    for sweep in range(steps) : 
        for i in range(vol) :
            k = randint(0, vol-1)
            h = S[nbhr(k)[0]] + S[nbhr(k)[1]] + S[nbhr(k)[2]] + S[nbhr(k)[3]]
            dE = 2*h*S[k]
            boltzmann = np.exp(-beta*dE)
            if ran(0., 1.) < boltzmann : S[k] = -S[k]
    return S 

@jit(nopython=True)
def calcMag(config) :
    ''' returns magentization for config '''
    return np.sum(config)

@jit(nopython=True)
def calcEnergy(config) :
    ''' returns the energy for config '''
    H = 0
    for i in range(vol) :
        H += config[i]*(config[nbhr(i)[0]] + config[nbhr(i)[1]] + config[nbhr(i)[2]] + config[nbhr(i)[3]])
    return H

file = open("data_L_modified_{}.txt".format(L), "w")
for t in range(len(betaJ)) :
    E1 = E2 = M1 = M2 = 0

    Ising(S, betaJ[t], 10000) # equilibrium steps chosen as 10000

    for i in range(MCsteps) : # count for calculating averages
        Ising(S, betaJ[t], 200)
        E1 += calcEnergy(S)
        E2 += calcEnergy(S)**2
        M1 += calcMag(S)
        M2 += calcMag(S)**2


    E[t] = n1*E1
    M[t] = n1*M1
    C[t] = (n1*E2 - n2*E1*E1)*(betaJ[t]**2)
    chi[t] = (n1*M2 - n2*M1*M1)*betaJ[t]
    file.write("% f\t% f\t% f\t% f\t% f\n" % (betaJ[t], E[t], M[t], C[t], chi[t]))
file.close()
        
f = plt.figure(figsize=(12, 9)) # plot the calculated values    

sp =  f.add_subplot(2, 2, 1 )
plt.scatter(betaJ, E, s=30, marker='o', c='r')
plt.xlabel(r"$\beta J$", fontsize=15)
plt.ylabel("Energy ", fontsize=15)
plt.axis('tight')
plt.grid(True)

sp =  f.add_subplot(2, 2, 2 )
plt.scatter(betaJ, abs(M), s=30, marker='o', c='b')
plt.xlabel(r"$\beta J$", fontsize=15)
plt.ylabel("Magnetization ", fontsize=15)
plt.axis('tight')
plt.grid(True)

sp =  f.add_subplot(2, 2, 3 )
plt.scatter(betaJ, C, s=30, marker='o', c='g')
plt.xlabel(r"$\beta J$", fontsize=15)
plt.ylabel("Specific Heat ", fontsize=15)
plt.axis('tight')
plt.grid(True)

sp =  f.add_subplot(2, 2, 4 )
plt.scatter(betaJ, chi, s=30, marker='o', c='orange')
plt.xlabel(r"$\beta J$", fontsize=15)
plt.ylabel("Susceptibility", fontsize=15)  
plt.axis('tight')
plt.grid(True)

plt.savefig("L_modified_{}".format(L))