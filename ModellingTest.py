import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

from ObservationModel import ObservationModel
import ObservationModelUtilities as omu

def modelColourExtraction(model):
    """ Returns list of colours for first two transmission curves
        passed to model after obervations have been run.
    Keyword arguments:
    model -- with relevant zRange,transDict,spectraDict    
    """
    
    if len(model.transCurves) < 2:
        raise ValueError("Model requires at least two defined Filters"
                         "for colour Calculation to be performed.")

    colourAtRedshiftList = []
    for key,mag in model.magnitudes.items():
        colour = mag[0] - mag[1]
        colourAtRedshiftList.append(colour)
    return colourAtRedshiftList        

def colourModelLyaResults(model,objectOverlayType="Default"):
    """ Running model for Lya of various EW values& plotting Results.
    Keyword arguments:
    model -- with relevant zRange,transDict,spectraDict
    objectOverlayType -- Default(For Lya),All,None.
    """
    equivalentWidthInuts = [50,100,150,200,250] #desired observable
    equivalentWidthInuts = [ew*2 for ew in equivalentWidthInuts]  
    colourDictOutputs = {} #key=EW

    for ewValue in equivalentWidthInuts:
        
        model.RestSpectrum.addEmissionLine(equivalentWidth=ewValue)
        model.makeObservedSpectra()
        model.runObservations()
        
        colourAtRedshiftList = modelColourExtraction(model)
        
        colourDictOutputs[ewValue] = colourAtRedshiftList
        model.RestSpectrum.removeEmissionLine('lya')

    ##loading in object data:
    ##file in format: ID,z,...,colour
    objectData = np.loadtxt("../potentialObjectsAllColours.ascii")

    ##For the y-Y colour vs z plot
    fig,axs = plt.subplots()
    colors = plt.cm.jet(np.linspace(0,1,len(equivalentWidthInuts)))
    for ewValue,color in zip(equivalentWidthInuts,colors):
        axs.plot(model.zRange,colourDictOutputs[ewValue],c=color,
                label="Obs. EW="+str(ewValue/2))
    axs.set_xlabel("z")
    axs.set_ylabel("y-Y Colour")
    axs.set_xlim(xmin=6,xmax=7.2)
    axs.set_ylim(ymin=-2,ymax=2)
    axs.axhline(0.0,color='k',ls=":")
    axs.legend(fontsize=8)

    omu.objectOverlayPlotter(axs,objectData,objectOverlayType)

    plt.show()

def colourModelIRACResults(model,objectOverlayType="None"):
    """ Running model for OIII line of various EW values& plotting Results.
    Keyword arguments:
    model -- with relevant zRange,transDict,spectraDict
    objectOverlayType -- Default(For Lya),All,None.
    """
    equivalentWidthInuts = [100,200,300,400,500] #desired observable
    colourDictOutputs = {} #key=EW

    for ewValue in equivalentWidthInuts:
        
        model.RestSpectrum.addEmissionLine(equivalentWidth=ewValue/3,
                                        lineCenter=4959,lineName='OIII-I')
        model.RestSpectrum.addEmissionLine(equivalentWidth=ewValue,
                                        lineCenter=5007,lineName="OIII-II")
        model.makeObservedSpectra()
        model.runObservations()
        
        colourAtRedshiftList = modelColourExtraction(model)
        
        colourDictOutputs[ewValue] = colourAtRedshiftList
        model.RestSpectrum.removeEmissionLine('OIII-I')
        model.RestSpectrum.removeEmissionLine('OIII-II')

    #loading in object data:
    #file in format: ID,z,y,Y,colour
    objectData = np.loadtxt("../derekYjhkMatchIracColours.ascii")

    fig,axs = plt.subplots()
    
    colors = plt.cm.jet(np.linspace(0,1,len(equivalentWidthInuts)))
    for ewValue,color in zip(equivalentWidthInuts,colors):
        axs.plot(model.zRange,colourDictOutputs[ewValue],c=color,
                label="{}$\AA$".format(ewValue))
    
    axs.set_xlabel("z")
    axs.set_ylabel(r"3.6$\mu$m - 4.5$\mu$m Colour")
    axs.set_xlim(xmin=5,xmax=9)
    axs.axhline(0.0,color='k',ls=":")
    legend = axs.legend(fontsize=8, title=r'EW for [OIII] 5007$\AA$ line:')
    plt.setp(legend.get_title(),fontsize=8)

    # To pass, say, every 2nd object, pass objectData[::2]
    print(objectData)
    omu.objectOverlayPlotter(axs,objectData,objectOverlayType)

    # ###Testing
    # iracObjectsFiles = ["../IRAC/RyanObjects_matchedToK_colour.ascii",
    #                     "../IRAC/RyanObjects_matchedToYJHK_colour.ascii"]
    # omu.iracObjectOverlayPlotter(axs,iracObjectsFiles)

    plt.show()


def main():

    zRange = [0,20]
    lineType = 'None'

    spectraDict = {"nuIndex":0.0,
                   "lowerLambda":100,
                   "upperLambda":30000,
                   "normalisationMag":25,
                   "activateIGM":True,
                   "activateLya":lineType}
    
    transDict = {"inputFile":["TransmissionCurveFiles/Subaru_HSC.Y.dat",
                              "TransmissionCurveFiles/Paranal_VISTA.Y.dat"]}

    
    model = ObservationModel(zRange,transDict,spectraDict)


    """Example usage & results of model."""
    # rossCalculations_v2 = np.loadtxt("ModellingDataFiles/Ross_delta_Y_colour_v2.dat")
    # gaussianCalulations = np.loadtxt("ModellingDataFiles/gaussian_delta_Y_colour.dat")
    # comparisonRoss = ("Ross Calc.",rossCalculations_v2,'k')
    # comparionsGaussian = ("Intr. EW=100",gaussianCalulations,'b')
    # omu.exampleModel(zRange,transDict,spectraDict,comparisonRoss,comparionsGaussian)

    """Run the animation of observational model."""
    # omu.runAnimation(model)

    """Plot results for the Lya similation w/ objects ontop plotted."""
    # colourModelLyaResults(model, objectOverlayType='All')


    """ Running the model for IRAC filters instead. """
    
    transDictIRAC = {"inputFile":["TransmissionCurveFiles/Spitzer_IRAC.I1.dat",
                                  "TransmissionCurveFiles/Spitzer_IRAC.I2.dat"]}
    modelIRAC = ObservationModel(zRange,transDictIRAC,spectraDict)
    colourModelIRACResults(modelIRAC,objectOverlayType="All")

if __name__ == '__main__':
    main()
