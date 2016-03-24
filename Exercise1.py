''' 
Multi-algorithm tutorial exercise #1 

Requires
- matplotlib
- numpy
- openpyxl

Recommended for interactive development
- ipython
- pprint
- pyreadline

Author: Jonathan Karr, karr@mssm.edu
Last updated: 3/17/2016
'''

#required libraries
from numpy import random
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as pyplot
import numpy as np
import utils #common code for exercises

#useful for pretty printing results (e.g. pp.print(...))
from pprint import PrettyPrinter
pp = PrettyPrinter(indent = 4)

'''Simulation parameters'''
modelFilename = 'Model.xlsx'
tMax = 10   #Length of simulation
tStep = 0.1 #Frequency at which to observe predicted cell state

'''Simulate cell starting from random initial conditions'''
def simulate(species, rxns):
    '''Initialize state'''
    #Time
    t = 0
    
    #Species copy numbers
    cnts = np.array(map(lambda x: np.round(random.normal(x['avg'], np.sqrt(x['avg']))), species))
	
	'''Initialize history'''
	tHist = np.reshape(np.linspace(0, tMax, num = tMax/tStep + 1), (1, -1))
	cntsHist = np.zeros((len(species), int(tMax/tStep + 1)))
	cntsHist[:, 0] = cnts
    
    '''Simulate dynamics'''
    #Hierarchical stochastic simulation algorithm
    while t < tMax:
        '''Calculate reactions and their propensities'''
        #Calculate reactions and their propensities for each individual submodel
        metRxns, metProps = calcMetabolismReactionsAndPropensities(cnts)
        trnRxns, trnProps = calcTranscriptionReactionsAndPropensities(cnts)
        trlRxns, trlProps = calcTranslationReactionsAndPropensities(cnts)
        cmpRxns, cmpProps = calcComplexationReactionsAndPropensities(cnts)
        dcyRxns, dcyProps = calcRnaDecayReactionsAndPropensities(cnts)

        #Calculate global reactions and propensities by concatenating reactions and 
        #propensities from individual submodels into single arrays
        rxns = np.concatenate(metRxns, trnRxns, trlRxns, cmpRxns, dcyRxns)
        props = np.concatenate(metProps, trnProps, trlProps, cmpProps, dcyProps)
        
        '''Select time to next reaction and the specific reaction'''
        #Select time to next reaction from exponential distribution
        dt = random.exponential(np.sum(props))
        
        #Select next reaction
        iRxn = random.multinomial(1, props)
        
        '''Update time and run selected reaction'''
        #Advance time
        t += dt
        
        #Run selected reaction
        specicesCnts += rxns[iRxn, :]
		
		'''Record state history'''
		cntsHist[:, np.ceil(t/tStep)] = cnts
	
	return (tHist, cntsHist)
        
'''Submodels'''
def calcMetabolismReactionsAndPropensities():
    return rxns, props
    
def calcTranscriptionReactionsAndPropensities():
    return rxns, props
    
def calcTranslationReactionsAndPropensities():
    return rxns, props
    
def calcComplexationReactionsAndPropensities(cnts):
    return rxns, props
    
def dcalcRnaDecayReactionsAndPropensities(cnts):
    return rxns, props
	
'''Plot'''
def plot(t, cnts, species, rxns):
	#reshape time array for plotting
	t = np.reshape(t, (-1, 1))
	
	#build dict of species indices
	speciesIdxs = {}
	for iSpecies in range(len(species)):
		speciesIdxs[species[iSpecies]['id']] = iSpecies

	'''NTPs'''
	#create figure
	fig = pyplot.figure()
	pyplot.plot(t, cnts[speciesIdxs['ATP'], :], label='ATP')
	pyplot.plot(t, cnts[speciesIdxs['CTP'], :], label='CTP')
	pyplot.plot(t, cnts[speciesIdxs['GTP'], :], label='GTP')
	pyplot.plot(t, cnts[speciesIdxs['UTP'], :], label='UTP')
	pyplot.title('NTPs')
	pyplot.xlabel('Time (s)')
	pyplot.ylabel('ATP')
	pyplot.legend()

	#save
	pp = PdfPages('ntps.pdf')
	pp.savefig(fig)
	pp.close()

'''Main function'''
if __name__ == "__main__":
	mdl = Model.getModelFromExcel(modelFilename)
    sim = Simulation.getSimulationFromModel(mdl)
    
    t, cnts = simulate(species, rxns)
	plot(t, cnts, species, rxns)