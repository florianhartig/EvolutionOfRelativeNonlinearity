'''
Created on 01.12.2014

@author: fberg
'''
from Invasion import *
import matplotlib.pyplot as plt
 
#Evolutionary and dynamic stability

#set Parameters
d = 0.05
r = 5
N1 = 200
N2 = 3
N = 100
K = 200
runs = 10
bvalues = np.logspace(-1, 1.3, 200)
bvalues.sort()
bposition = len(bvalues)
invasion = Invasion()
 
# Create result Arrays
result = np.zeros((200, 200))

# Run invasion  
for a in range(0, bposition):
   b1 = bvalues[a]
   for b in range(0, bposition):
       b2 = bvalues[b]
       result[a, b] = invasion.runInvasion(10, 500, 500, N1, N2, r, d, K, b1, b2)
        

# Create Plot 'Invasion probability'
y, x = np.meshgrid(bvalues, bvalues)
z = result

plt.pcolor(x, y, z, cmap='Greys')
plt.colorbar()
plt.xscale('log')
plt.yscale('log')
plt.xlim(0, 19.95)
plt.ylim(0, 19.95)
plt.axvline(x=2.533333333, color='k', linestyle='--')
plt.axhline(y=2.533333333, color='k', linestyle='--')
tick = [0.5, 1, 2, 3, 4, 5, 6, 7, 8, 10, 15]
plt.xticks(tick, tick)
plt.yticks(tick, tick)
plt.xlabel('Resident density-compensation strategy b')
plt.ylabel('Invading density-compensation strategy b')
plt.title('(a) Invasion probability')
plt.savefig('Result3a.png')

# Clear Result arrays
result = np.zeros((200,200))
result2 = np.zeros((200,200))

# Run coexistence
for a in range(0, bposition):
   b1 = bvalues[a]
   for b in range(0, bposition):
       b2 = bvalues[b]
       resulta=invasion.runCoexistence(N, r, d, K, b1, b2)
       result[a,b] = resulta[0]
       result2[a,b] = resulta[1]


# Create Graphic 'Log10 time until competitive exclusion'
y, x = np.meshgrid(bvalues, bvalues)
z = result

plt.clf()
plt.pcolor(x, y, z, cmap='jet')
plt.colorbar()
plt.xscale('log')
plt.yscale('log')
plt.xlim(0,19.95)
plt.ylim(0,19.95)
plt.xticks(tick,tick)
plt.yticks(tick,tick)
plt.xlabel('Resident density-compensation strategy b')
plt.ylabel('Invading density-compensation strategy b')
plt.title('(b) Log10 time until competitive exclusion')
plt.savefig('Result4a.png')
  
# Create Graphic 'Surviving species of competitive exclusion'
z = result2
  
plt.clf()
plt.pcolor(x, y, z, cmap='Greys')
plt.xscale('log')
plt.yscale('log')
plt.xlim(0, 19.95)
plt.ylim(0, 19.95)
plt.xticks(tick, tick)
plt.yticks(tick, tick)
plt.xlabel('Resident density-compensation strategy b')
plt.ylabel('Invading density-compensation strategy b')
plt.title('(c) Surviving species of competitive exclusion')
plt.savefig('Result4a2.png')

# Development of coexisting communities
# Set Parameters
N = 500
K = 1000
m = 0.3
b1 = 0.9
b2 = 5.4
 
# Run Mutationdevelopment for coexisting community
bValues = np.zeros((1500,3))
bValues = invasion.runMutationdevelopment_coexistingCommunity(N, r, d, K, b1, b2, m)
  
# Create Graphic 'Collapse of a coexisting community'
plt.clf()
plt.plot(bValues[:,0],bValues[:,1],marker='.',linestyle='')
plt.plot(bValues[:,0],bValues[:,2],marker='.',linestyle='',color='r')
plt.xlabel('Time')
plt.ylabel('Density-compensation strategy b')
plt.title('(a) Collapse of a coexisting community')
plt.axvline(x=30000, color='k', linestyle='--')
plt.savefig('Result5.png')

# Change of b-values
b1 = 0.5
b2 = 6
 
# Run Mutationdevelopment for independent Populations
bValues = invasion.runMutationdevelopment_independentPopulations(N, r, d, K, b1, b2, m)
  
# Create Graphic 'Convergent evolution of independent populations'
plt.clf()
plt.plot(bValues[:,0],bValues[:,1],marker='.',linestyle='')
plt.plot(bValues[:,0],bValues[:,2],marker='.',linestyle='',color='r')
plt.xlabel('Time')
plt.ylabel('Density-compensation strategy b')
plt.title('(b) Convergent evolution of independent populations')
plt.axvline(x=30000, color='k', linestyle='--')
plt.savefig('Result6.png')

# ESS
# Set Parameters
r_values = np.arange(1,6,0.0125)
r_position = 0
d_values = np.arange(0,0.5,0.00125)
d_position = 0
N = 200
K = 200
m = 0.1

# Create result Arrays
result1 = np.zeros((400,400))
result2 = np.zeros((400,400))

# Run Population development
for r in r_values:
    d_position = 0
    for d in d_values:
        pop = Population(m,mutationON=True)
        pop.addNdifferentIndividuals(N, r, d, K, resident=False)
        pop.updatePopulation(10000)
        result1[r_position,d_position]=pop.getAverageB()
        result2[r_position,d_position]=pop.getCriticalB(r, d)-pop.getAverageB()
        d_position += 1
    r_position += 1
    
#Create Graphic 'ESS of individual-based model'
y, x = np.meshgrid(d_values,r_values)
z = result1

plt.clf()
plt.pcolor(x, y, z, cmap='hot', vmin=0, vmax=15)
plt.colorbar()
plt.xlabel('Intrinsic growth rate r')
plt.ylabel('Density-independent mortality d')
plt.title('(a) ESS of individual-based model')
plt.savefig("Result1.png")

# Create Graphic 'Deviation of ESS from critical b-value'
z = result2

plt.clf()
plt.pcolor(x,y,z,cmap='seismic',vmin=-10,vmax=10)
plt.colorbar()
plt.xlabel('Intrinsic growth rate r')
plt.ylabel('Density-independent mortality d')
plt.title('(b) Deviation of ESS from critical b-value')
plt.savefig("Result2.png")
