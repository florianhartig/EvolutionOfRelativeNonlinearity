'''
Created on 18.12.2014

@author: fberg
'''
from Population import *
import math as ma

class Invasion:
    
    def __init__(self):
        None
          
    def runInvasion(self, repetitions, delay, dataTime, N1, N2, r, d, K, b1, b2):
    #: Function for running the invasion 
    #: @param repetitions: repetitions
    #: @return: returns the fraction of runs  in which  the invader was successful
        invasion = 0
        for rep in range(0,repetitions):
            population = Population(0, mutationON=False)
            population.addIndividual(r, d, N1, b1, K, resident=True)
            population.updatePopulation(delay)
            population.addIndividual(r, d, N2, b2, K, resident=False)
            population.updatePopulation(dataTime)
            if population.invaderPresent(): invasion += 1.0
        return invasion / repetitions

    def runCoexistence(self, N, r, d, K, b1, b2):
        timeUntilDeath = 0 
        survivor = 0
        population = Population(0, mutationON=False)
        population.addIndividual(r, d, N, b1, K, resident=True)
        population.addIndividual(r, d, N, b2, K, resident=False)
        while population.residentinvaderPresent():
            population.updatePopulation(1)
            timeUntilDeath += 1
        if population.residentPresent(): survivor = 1
        if population.invaderPresent(): survivor = 2
        return ma.log10(timeUntilDeath), survivor
    
    def runMutationdevelopment_coexistingCommunity(self, N, r, d, K, b1, b2, m):
        bValues = np.zeros((15000,3))
        time = 1
        resultPosition = 0
        population = Population(m, mutationON=False)
        population.addIndividual(r, d, N, b1, K, resident=True)
        population.addIndividual(r, d, N, b2, K, resident=False)
        while time <= 30000:
            population.updatePopulation(1)
            if time % 100 == 1: 
                bValues[resultPosition]= [time, population.getResidentAverageB(),population.getInvaderAverageB()]
                resultPosition += 1
            time += 1
        population.mutationON = True
        while time <= 150000:
            population.updatePopulation(1)
            if time % 100 == 1: 
                bValues[resultPosition] = [time, population.getResidentAverageB(),population.getInvaderAverageB()]
                resultPosition += 1
            time += 1
        return bValues
    
    def runMutationdevelopment_independentPopulations(self, N, r, d, K, b1, b2, m):
        bValues = np.zeros((15000,3))
        time = 1
        resultPosition = 0
        population1 = Population(m, mutationON=False)
        population1.addIndividual(r, d, N, b1, K, resident=True)
        population2 = Population(m, mutationON=False)
        population2.addIndividual(r, d, N, b2, K, resident=True)
        while time <= 30000:
            population1.updatePopulation(1)
            population2.updatePopulation(1)
            if time % 100 == 1: 
                bValues[resultPosition] = [time, population1.getAverageB(), population2.getAverageB()]
                resultPosition += 1
            time += 1
        population1.mutationON = True
        population2.mutationON = True
        while time <= 150000:
            population1.updatePopulation(1)
            population2.updatePopulation(1)
            if time % 100 == 1: 
                bValues[resultPosition] = [time, population1.getAverageB(), population2.getAverageB()]
                resultPosition += 1
            time += 1
        return bValues