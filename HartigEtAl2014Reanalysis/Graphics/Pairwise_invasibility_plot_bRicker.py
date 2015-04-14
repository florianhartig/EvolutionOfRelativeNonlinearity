'''
Created on 05.12.2014

@author: fberg
'''
import numpy as np
import math as ma
import matplotlib.pyplot as plt
from matplotlib.pyplot import savefig

# Definition of the Model
from PopulationModels import Ricker

#Creates resident population
def residentPopulation (N, cycles, d, r, K, b):
    resident = []
    T = cycles +1
    for t in range (0,T):
        rep = Ricker.getReproductionRatio(r, d, N, b, K)
        N = N*rep
        if t >= 100:
            resident += [N,]
    return resident

#Calculates fitness of the Invader
def invaderFitness (resident, d, r, K, b):
    fit = 0
    T = len(resident)
    for t in range (0,T):
        N = resident[t]
        rep = Ricker.getReproductionRatio(r, d, N, b, K)
        if rep != 0: fit += ma.log(rep)
    fitness = fit/T
    return fitness

#Sets initial Values
d = 0.05
r = 0.5
N = 100
K = 200
cycles = 500
bvalues = np.logspace(-1,1.3,200)
bvalues.sort()
bposition = len(bvalues)

#Create result Arrays
result = np.zeros((200,200))
result2 = np.zeros((200,200))

#Fill result Array
for a in range (0,bposition):
    bR = bvalues[a]
    resident = residentPopulation(N, cycles, d, r, K, bR)
    for b in range(0,bposition):
        bI = bvalues[b]
        fitness = invaderFitness(resident, d, r, K, bI)
        result[a,b] = fitness
        
#Fill result Array 2
for a in range (0,bposition):
    for b in range (0,bposition):
        if (result[a,b] > 0) and (result[b,a] > 0):
            result2[a,b] = 1
            result2[b,a] = 1
        else:
            result2[a,b] = 0
            result2[b,a] = 0
    
#Creates plot data
y, x = np.meshgrid(bvalues,bvalues)
z = result
z2 = result2

#Define colormap
cdict3 = {'red':  ((0.0, 0.0, 0.0),
                   (0.25,0.0, 0.0),
                   (0.5, 0.8, 1.0),
                   (0.75,1.0, 1.0),
                   (1.0, 0.4, 1.0)),

         'green': ((0.0, 0.0, 0.0),
                   (0.25,0.0, 0.0),
                   (0.5, 0.9, 0.9),
                   (0.75,0.0, 0.0),
                   (1.0, 0.0, 0.0)),

         'blue':  ((0.0, 0.0, 0.4),
                   (0.25,1.0, 1.0),
                   (0.5, 1.0, 0.8),
                   (0.75,0.0, 0.0),
                   (1.0, 0.0, 0.0))
        }
plt.register_cmap(name='BlueRed', data=cdict3)

#Defines the graphic 1
plt.pcolor(x,y,z,cmap='BlueRed', vmin = -2, vmax= 2)
plt.colorbar()
plt.xscale('log')
plt.yscale('log')
plt.xlim(0,19.95)
plt.ylim(0,19.95)
plt.axvline(x=4.45725, color='k', linestyle='--')
tick = [0.5,1,2,3,4,5,6,7,8,10,15]
plt.xticks(tick,tick)
plt.yticks(tick,tick)
plt.title('(b) Ricker model ')
plt.xlabel('Resident density-compensation strategy b')
plt.ylabel('Invading density-compensation strategy b')
savefig("PIP_Ricker1.png")

#Defines graphic 2
plt.clf()
plt.pcolor(x,y,z2, cmap= "Greys")
plt.xscale('log')
plt.yscale('log')
plt.xlim(0,19.95)
plt.ylim(0,19.95)
plt.xticks(tick,tick)
plt.yticks(tick,tick)
plt.title('(b) Ricker model ')
plt.xlabel('Resident density-compensation strategy b')
plt.ylabel('Invading density-compensation strategy b')
savefig("PIP_Ricker2.png")

print "done"
