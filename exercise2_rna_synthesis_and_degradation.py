#!/usr/bin/python

'''
Simulates transcription and RNA degradation submodels

@author Jonathan Karr, karr@mssm.edu
@date 3/24/2016
'''

#required libraries
from model import getModelFromExcel, SsaSubmodel #code for model in exercises
from numpy import random
import analysis #code to analyze simulation results in exercises
import numpy as np
import os
import util

#simulation parameters
MODEL_FILENAME = 'data/Model.xlsx'
TIME_STEP = 10 #time step on simulation (s)
TIME_STEP_RECORD = TIME_STEP #Frequency at which to observe predicted cell state (s)
OUTPUT_DIRECTORY = 'out/exercise2_rna_synthesis_and_degradation'
RANDOM_SEED = 0

#simulates model
def simulate(model):    
    #Get transcription, RNA degradation submodels
    transcriptionSubmodel = model.getComponentById('Transcription')
    rnaDegradationSubmodel = model.getComponentById('RnaDegradation')

    #get parameters
    cellCycleLength = model.getComponentById('cellCycleLength').value
    
    #seed random number generator to generate reproducible results
    random.seed(RANDOM_SEED)

    #Initialize state
    model.calcInitialConditions()

    time = 0 #(s)
    volume = model.volume
    extracellularVolume = model.extracellularVolume
    speciesCounts = model.speciesCounts

    #get data to mock other submodels
    metabolismSubmodel = model.getComponentById('Metabolism')   
    netMetabolismReaction = np.zeros((len(model.species), len(model.compartments)))
    for part in metabolismSubmodel.getComponentById('MetabolismProduction').participants:
        netMetabolismReaction[part.species.index, part.compartment.index] = -part.coefficient
            
    translationSubmodel = model.getComponentById('Translation')   
    netTranslationReaction = np.zeros((len(model.species), len(model.compartments)))
    for rxn in translationSubmodel.reactions:
        for part in rxn.participants:
            if part.species.type == 'Protein':
                initCopyNumber = model.speciesCounts[part.species.index, part.compartment.index]
        for part in rxn.participants:
            netTranslationReaction[part.species.index, part.compartment.index] += part.coefficient * initCopyNumber

    #Initialize history
    timeMax = cellCycleLength #(s)
    nTimeSteps = int(timeMax / TIME_STEP + 1)
    nTimeStepsRecord = int(timeMax / TIME_STEP_RECORD + 1)
    timeHist = np.linspace(0, timeMax, num = nTimeStepsRecord)

    volumeHist = np.full(nTimeStepsRecord, np.nan)
    volumeHist[0] = volume

    extracellularVolumeHist = np.full(nTimeStepsRecord, np.nan)
    extracellularVolumeHist[0] = extracellularVolume

    speciesCountsHist = np.zeros((len(model.species), len(model.compartments), nTimeStepsRecord))
    speciesCountsHist[:, :, 0] = speciesCounts
            
    #Simulate dynamics
    print 'Simulating for %d time steps from 0-%d s' % (nTimeSteps, timeMax)
    for iTime in range(1, nTimeSteps):
        time = iTime * TIME_STEP
        if iTime % 100 == 1:
            print '\tStep = %d, t=%.1f s' % (iTime, time)
        
        #simulate submodels
        model.setSpeciesCountsDict(SsaSubmodel.stochasticSimulationAlgorithm(
            model.getSpeciesCountsDict(), 
            model.getSpeciesVolumesDict(), 
            [transcriptionSubmodel.reactions, rnaDegradationSubmodel.reactions], 
            model.volume,
            TIME_STEP,
            ))
        
        #mock other submodels
        model.speciesCounts += netMetabolismReaction  * np.log(2) / cellCycleLength * np.exp(np.log(2) * time / cellCycleLength) * TIME_STEP
        model.speciesCounts += netTranslationReaction * np.log(2) / cellCycleLength * np.exp(np.log(2) * time / cellCycleLength) * TIME_STEP
        
        #update mass, volume        
        model.calcMass()
        model.calcVolume()
                
        #Record state
        volumeHist[iTime] = model.volume
        extracellularVolumeHist[iTime] = model.extracellularVolume
        speciesCountsHist[:, :, iTime] = model.speciesCounts
    
    return (timeHist, volumeHist, extracellularVolumeHist, speciesCountsHist)
    
#plot results
def analyzeResults(model, time, volume, extracellularVolume, speciesCounts):
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
        speciesCounts = speciesCounts, 
        units = 'molecules',
        selectedSpeciesCompartments = ['RnaPolymerase-Rna[c]', 'Adk-Rna[c]', 'Apt-Rna[c]', 'Cmk-Rna[c]'], 
        fileName = os.path.join(OUTPUT_DIRECTORY, 'RNA.pdf')
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
    time, volume, extracellularVolume, speciesCounts = simulate(model)
    analyzeResults(model, time, volume, extracellularVolume, speciesCounts)
