'''
Created on 27.11.2014

@author: Franziska Berg
'''
from __future__ import division

import numpy as np
from copy import deepcopy
from PopulationModels import MSS

class Individual:
    
    def __init__(self,r,d,N,b,K,resident=True):
        self.number = N
        self.r = r
        self.d = d
        self.K = K
        self.resident = resident
        self.b_value = b
       
    def reproduce(self,N):
        poisson = self.number*MSS.getReproductionRatio(self.r, self.d, N, self.b_value, self.K)    
        #self.r*(1-self.d)/(1+(self.r-1)*(N/self.K)**self.b_value)
        self.number = np.random.poisson(poisson)
        
   
    def getMutant(self):
        self.number -= 1
        newMutant = deepcopy(self)
        newMutant.number = 1
        newMutant.b_value = self.b_value + np.random.normal(scale=0.15)
        if newMutant.b_value < 0.17:
            newMutant.b_value = 0.17
        return(newMutant)
        
    def getIndividualNumber(self):
        return self.number
    
    def getB(self):
        return self.b_value
    
            