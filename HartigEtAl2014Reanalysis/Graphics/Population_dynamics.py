'''
Created on 11.11.2014

@author: fberg
'''
import matplotlib.pyplot as plt
import numpy as np

r = 5
d = 0.05
N = 100
K = 200

from PopulationModels import MSS

result = np.zeros((280*400,3),"float")
s = 0

for b in np.linspace(1,7,280):
    for i in range (1,500):
        Nt = N*MSS.getReproductionRatio(r, d, N, b, K)
        if i > 100:
            result[s,:] =([b, N, K])
            s+= 1
        N = Nt
       
plt.plot(result[:,0],result[:,1]/result[:,2],'k.')
plt.xlabel('Density- compensation strategy b')
plt.ylabel('Population size / carrying capacity N/K')
plt.axis([1,7,0.,2.5])
plt.axvline(x=2.533333333, color='k', linestyle='--')
plt.show()