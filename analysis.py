''' 
Analysis utility functions

@author Jonathan Karr, karr@mssm.edu
@date 3/26/2016
'''

#required libraries
from matplotlib.backends.backend_pdf import PdfPages
from model import Submodel
from util import N_AVOGADRO
import matplotlib.pyplot as pyplot
import numpy as np
import re

def plot(model, time = np.zeros(0), 
    speciesCounts = None, volume = np.zeros(0), extracellularVolume = np.zeros(0),
    selectedSpecies = [], units = 'mM', title = '', fileName = ''):

    #reshape time array for plotting and convert to hours
    time = np.reshape(time.copy(), (-1, 1)) / 3600
    
    #create figure
    fig = pyplot.figure()

    #plot results
    yMin = 1e12
    yMax = -1e12
    for speciesId in selectedSpecies:
        #extract data
        match = re.match('^(?P<speciesId>[a-z0-9\-_]+)\[(?P<compartmentId>[a-z0-9\-_]+)\]$', speciesId, re.I).groupdict()
        compartmentId = match['compartmentId']

        yData = speciesCounts[speciesId]
            
        #scale
        if compartmentId == 'c':
            V = volume
        else:
            V = extracellularVolume
        
        if units == 'pM':
            scale = 1 / N_AVOGADRO / V * 1e12
        elif units == 'nM':
            scale = 1 / N_AVOGADRO / V * 1e9
        elif units == 'uM':
            scale = 1 / N_AVOGADRO / V * 1e6
        elif units == 'mM':
            scale = 1 / N_AVOGADRO / V * 1e3
        elif units == 'M':
            scale = 1 / N_AVOGADRO / V * 1e0
        elif units == 'molecules':
            scale = 1
        else:
            raise Exception('Invalid units "%s"' % units)
            
        yData *= scale
        
        #update range
        yMin = min(yMin, np.min(yData))
        yMax = max(yMax, np.max(yData))

        #add to plot
        pyplot.plot(time, yData, label=speciesId)
        
    #set axis limits
    pyplot.xlim((0, time[-1]))
    pyplot.ylim((yMin, yMax))
    
    #add axis labels and legend
    if title:
        pyplot.title(title)
    
    pyplot.xlabel('Time (h)')

    if units == 'molecules':
        pyplot.ylabel('Copy number')
    else:
        pyplot.ylabel('Concentration (%s)' % units)
    
    if len(selectedSpecies) > 1:
        pyplot.legend()

    #save
    if fileName:
        pp = PdfPages(fileName)
        pp.savefig(fig)
        pp.close()
