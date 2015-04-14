'''
Created on 22.01.2015

@author: fberg
'''
from __future__ import division
import math as ma

class PopulationModel(object):
    '''
    classdocs
    '''

    def __init__(self, params):
        '''
        Constructor
        '''
    @staticmethod  
    def getReproductionRatio():
        raise NotImplementedError
    
    @staticmethod
    def plotFunctionalResponse():
        raise NotImplementedError
    
    
class MSS(PopulationModel):
    '''
    classdocs
    '''

    def __init__(self, params):
        '''
        Constructor
        '''
    @staticmethod      
    def getReproductionRatio(r,d,N,b,K):
        rep = r*(1-d)/(1+(r-1)*(N/K)**b)
        return rep
    
    @staticmethod
    def getReproductionRatio2(r,d,pop,b):
        rep = r*(1-d)/(1+(r-1)*(pop)**b)
        return rep
        
            
    
class Ricker(PopulationModel):
    
    def __init__(self):
        '''
        '''
        
    @staticmethod
    def getReproductionRatio(r,d,N,b,K):
        rep = (1-d)*ma.exp(r*(1-(N/K)**b))
        return rep
    
    
class Hassel(PopulationModel):
    
    def __init__(self):
        '''
        '''
        
    @staticmethod
    def getReproductionRatio(r, d, N, b, K):
        rep = (1-d)*r/(1+(r**(1/b)-1)*N/K)**b
        return rep
