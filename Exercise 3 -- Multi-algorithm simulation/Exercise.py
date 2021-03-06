#!/usr/bin/python

'''
Simulates metabolism submodel

@author Jonathan Karr, karr@mssm.edu
@date 3/24/2016
'''

#required libraries
from model import getModelFromExcel, Submodel, SsaSubmodel #code for model in exercises
from numpy import random
from util import N_AVOGADRO
import analysis #code to analyze simulation results in exercises
import numpy as np
import os

#simulation parameters
MODEL_FILENAME = 'Model.xlsx'
TIME_STEP = 10 #time step on simulation (s)
TIME_STEP_RECORD = TIME_STEP #Frequency at which to observe predicted cell state (s)
OUTPUT_DIRECTORY = '.'
RANDOM_SEED = 10000000

#simulates model
def simulate(model):
    #Get FBA, SSA submodels
    ssaSubmodels = []
    for submodel in model.submodels:
        if isinstance(submodel, SsaSubmodel):
            ssaSubmodels.append(submodel)
            
    metabolismSubmodel = model.getComponentById('Metabolism')

    #get parameters
    cellCycleLength = model.getComponentById('cellCycleLength').value

    #seed random number generator to generate reproducible results
    random.seed(RANDOM_SEED)

    #Initialize state
    model.calcInitialConditions()

    time = 0 #(s)

    #Initialize history
    timeMax = cellCycleLength #(s)
    nTimeSteps = int(timeMax / TIME_STEP + 1)
    nTimeStepsRecord = int(timeMax / TIME_STEP_RECORD + 1)
    timeHist = np.linspace(0, timeMax, num = nTimeStepsRecord)

    volumeHist = np.full(nTimeStepsRecord, np.nan)
    volumeHist[0] = model.volume

    growthHist = np.full(nTimeStepsRecord, np.nan)
    growthHist[0] = np.log(2) / cellCycleLength

    speciesCountsHist = np.zeros((len(model.species), len(model.compartments), nTimeStepsRecord))
    speciesCountsHist[:, :, 0] = model.speciesCounts
            
    #Simulate dynamics
    print 'Simulating for %d time steps from 0-%d s' % (nTimeSteps, timeMax)
    for iTime in range(1, nTimeSteps):
        time = iTime * TIME_STEP
        if iTime % 100 == 1:
            print '\tStep = %d, t = %.1f s' % (iTime, time)
        
        #simulate submodels
        #<Synchronize metabolism submodel: submodel.updateLocalCellState(model)>
        #<Calculate reaction bounds: submodel.calcReactionBounds(TIME_STEP)>
        #<Optimize FBA: submodel.calcReactionFluxes(TIME_STEP)>
        #<Update metabolites: submodel.updateMetabolites(TIME_STEP)>
        #<Synchronize metabolism submodel: submodel.updateGlobalCellState(model)>

        speciesCountsDict = model.getSpeciesCountsDict()        
        time2 = 0
        while time2 < TIME_STEP:
            time = 0
            
            #calculate concentrations
            #<convert species counts to concentrations>
        
            #calculate propensities
            #<iterate over SSA submodels and calculate reaction rates: Submodel.calcReactionRates(submodel.reactions, specicesConcentrations)>
            
            #Select time to next reaction from exponential distribution
            #<dt = random.exponential(...)>
            
            #Select next reaction
            #<iSubmodel = random.choice(...)>
            #<iRxn = random.choice(...)>

            #update time
            #<time2 += dt>
            
            #execute reaction
            #<selectedSubmodel = ...>
            #<selectedReaction = ...>
            #<selectedSubmodel.executeReaction(selectedReaction)>

        model.setSpeciesCountsDict(speciesCountsDict)
        
        #update mass, volume        
        model.calcMass()
        model.calcVolume()
                
        #Record state
        volumeHist[iTime] = model.volume
        growthHist[iTime] = model.growth
        speciesCountsHist[:, :, iTime] = model.speciesCounts

    return (timeHist, volumeHist, growthHist, speciesCountsHist)
    
#plot results
def analyzeResults(model, time, volume, growth, speciesCounts):
    if not os.path.exists(OUTPUT_DIRECTORY):
        os.makedirs(OUTPUT_DIRECTORY)
        
    cellComp = model.getComponentById('c')
    
    totalRna = np.zeros(len(time))
    totalProt = np.zeros(len(time))
    for species in model.species:
        if species.type == 'RNA':
            totalRna += speciesCounts[species.index, cellComp.index, :]
        elif species.type == 'Protein':
            totalProt += speciesCounts[species.index, cellComp.index, :]
	
    analysis.plot(
        model = model, 
        time = time, 
        yDatas = {'Volume': volume},
        fileName = os.path.join(OUTPUT_DIRECTORY, 'Volume.pdf')
        )
        
    analysis.plot(
        model = model, 
        time = time, 
        yDatas = {'Growth': growth},
        fileName = os.path.join(OUTPUT_DIRECTORY, 'Growth.pdf')
        )
            
    analysis.plot(
        model = model, 
        time = time, 
        yDatas = {'RNA': totalRna},
        fileName = os.path.join(OUTPUT_DIRECTORY, 'Total RNA.pdf')
        )
        
    analysis.plot(
        model = model, 
        time = time, 
        yDatas = {'Protein': totalProt},
        fileName = os.path.join(OUTPUT_DIRECTORY, 'Total protein.pdf')
        )
    
    analysis.plot(
        model = model, 
        time = time, 
        volume = volume, 
        speciesCounts = speciesCounts, 
        units = 'molecules',
        selectedSpeciesCompartments = ['ATP[c]', 'CTP[c]', 'GTP[c]', 'UTP[c]'], 
        fileName = os.path.join(OUTPUT_DIRECTORY, 'NTPs.pdf')
        )

    analysis.plot(
        model = model, 
        time = time, 
        volume = volume, 
        speciesCounts = speciesCounts, 
        selectedSpeciesCompartments = ['AMP[c]', 'CMP[c]', 'GMP[c]', 'UMP[c]'], 
        units = 'uM',
        fileName = os.path.join(OUTPUT_DIRECTORY, 'NMPs.pdf')
        )
        
    analysis.plot(
        model = model, 
        time = time, 
        volume = volume,
        speciesCounts = speciesCounts, 
        selectedSpeciesCompartments = ['ALA[c]', 'ARG[c]', 'ASN[c]', 'ASP[c]'], 
        units = 'uM',
        fileName = os.path.join(OUTPUT_DIRECTORY, 'Amino acids.pdf')
        )
        
    analysis.plot(
        model = model, 
        time = time, 
        speciesCounts = speciesCounts, 
        units = 'molecules',
        selectedSpeciesCompartments = ['RnaPolymerase-Protein[c]', 'Adk-Protein[c]', 'Apt-Protein[c]', 'Cmk-Protein[c]'], 
        fileName = os.path.join(OUTPUT_DIRECTORY, 'Proteins.pdf')
        )
        
#main function
if __name__ == "__main__":
    model = getModelFromExcel(MODEL_FILENAME)
    time, volume, growth, speciesCounts = simulate(model)
    analyzeResults(model, time, volume, growth, speciesCounts)
    
    #Check if simulation implemented correctly
    volumeChange = (volume[-1] - volume[0]) / volume[0]
    if volumeChange < 0.9 or volumeChange > 1.1:
        raise Exception('Volume should approximately double over the simulation.')
