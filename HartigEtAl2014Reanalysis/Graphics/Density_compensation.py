'''
Created on 11.11.2014

@author: fberg
'''
import numpy as np
import matplotlib.pyplot as plt
from PopulationModels import MSS

d = 0.05
r = 5
blist = [0.1,1,2.5,6]

for b in blist:
    popList=[]
    repList=[]
    
    for pop in np.arange(0.,2.,0.1):
        popList.append(pop)
        rep = MSS.getReproductionRatio2(r, d, pop, b)
        repList.append(rep)
    
    plt.plot(popList,repList, label = 'b =' + str(b))
    plt.legend()
    
plt.xlabel('Population size / carrying capacity N/K')
plt.ylabel('Reproduction ratio N(t+1)/N(t)')
plt.axhline(y=1, color='k', linestyle='--')
plt.axvline(x=1, color='k', linestyle='--')
plt.show()
