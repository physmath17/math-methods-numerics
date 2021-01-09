import numpy as np
from numpy.random import rand
import matplotlib.pyplot as plt
from numba import jit
from random import uniform as ran, choice

class Ising():
    ''' Simulating the Ising model '''    
    L, temp = 10, .4        # Initialse the lattice
    vol = L*L
    config = np.array([choice([-1, 1]) for k in range(vol)])
    ## monte carlo moves
    @jit(nopython=True)
    def MCmove(self, beta):
        ''' This is to execute the Monte Carlo moves using 
        Metropolis algorithm such that detailed
        balance condition is satisified'''
        for k in range(vol) :
            h = sum(config[self.nbhr(k)[j]] for j in range(4))
            dE = 2*h*config[k]
            boltzmann = np.exp(-beta*dE)
            if ran(0., 1.) < boltzmann : config[k] = -config[k]
        return config
    
    @jit(nopython=True)
    def simulate(self):   
        ''' This module simulates the Ising model'''
        f = plt.figure(figsize=(15, 15), dpi=80);    
        self.configPlot(f, config, 0, vol, 1);
        
        msrmnt = 1001
        for i in range(msrmnt):
            self.MCmove(config, vol, 1.0/temp)
            if i == 1 :       self.configPlot(f, config, i, vol, 2);
            if i == 4 :       self.configPlot(f, config, i, vol, 3);
            if i == 50 :      self.configPlot(f, config, i, vol, 4);
            if i == 100 :     self.configPlot(f, config, i, vol, 5);
            if i == 200 :     self.configPlot(f, config, i, vol, 5);
            if i == 500 :     self.configPlot(f, config, i, vol, 5);
            if i == 1000 :    self.configPlot(f, config, i, vol, 6);

    @jit(nopython=True)
    def nbhr(self, k) :
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

    def configPlot(self, f, config, i, vol, n_):
        ''' This modules plts the configuration once passed to it along with time etc '''
        X, Y = np.meshgrid(range(N), range(N))
        sp =  f.add_subplot(3, 3, n_ )  
        plt.setp(sp.get_yticklabels(), visible=False)
        plt.setp(sp.get_xticklabels(), visible=False)      
        plt.pcolormesh(X, Y, config, cmap=plt.cm.RdBu);
        plt.title('Time=%d'%i); plt.axis('tight')    
    plt.show()

kj = Ising()
kj.simulate()