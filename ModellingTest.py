import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

from ObservationModel import ObservationModel
import ObservationModelUtilities as omu


def colourModelLyaResults(model,objectOverlayType="Default"):
    """Running model for Lya of various EW values& plotting Results.
    Keyword arguments:
    model -- with relevant zRange,transDict,spectraDict
    objectOverlayType -- Default,All,None.
    """
    equivalentWidthInuts = [50,100,150,200,250] #desired observable
    equivalentWidthInuts = [ew*2 for ew in equivalentWidthInuts]  
    colourDictOutputs = {} #key=EW

    for ewValue in equivalentWidthInuts:
        
        model.RestSpectrum.addEmissionLine(equivalentWidth=ewValue)
        model.makeObservedSpectra()
        model.runObservations()
        
        colourAtRedshiftList = []
        
        for key,mag in model.magnitudes.items():
            yminusY = mag[0] - mag[1]
            colourAtRedshiftList.append(yminusY)
        
        colourDictOutputs[ewValue] = colourAtRedshiftList
        model.RestSpectrum.removeEmissionLine('lya')

    ##loading in object data:
    ##file in format: ID,z,y,Y,colour
    objectData = np.loadtxt("../potentialObjectsColours.ascii")

    ##For the y-Y colour vs z plot
    fig,axs = plt.subplots()
    colors = plt.cm.jet(np.linspace(0,1,len(equivalentWidthInuts)))
    for ewValue,color in zip(equivalentWidthInuts,colors):
        axs.plot(model.zRange,colourDictOutputs[ewValue],c=color,
                label="Obs. EW="+str(ewValue/2))
    axs.set_xlabel("z")
    axs.set_ylabel("y-Y Colour")
    axs.set_xlim(xmin=6,xmax=7.5)
    axs.set_ylim(ymin=-2,ymax=2)
    axs.axhline(0.0,color='k',ls=":")
    axs.legend(fontsize=8)

    omu.objectOverlayPlotter(axs,objectData,objectOverlayType)

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
    # rossCalculations_v2 = np.loadtxt("Ross_delta_Y_colour_v2.dat")
    # gaussianCalulations = np.loadtxt("gaussian_delta_Y_colour.dat")
    # comparisonRoss = ("Ross Calc.",rossCalculations_v2,'k')
    # comparionsGaussian = ("Intr. EW=100",gaussianCalulations,'b')
    # omu.exampleModel(zRange,transDict,spectraDict,comparisonRoss,comparionsGaussian)

    """Run the animation of observational model."""
    # omu.runAnimation(model)

    """Plot results for the Lya similation w/ objects ontop plotted."""
    # colourModelLyaResults(model)


    """ Running the model for IRAC filters instead. """
    # transDictIRAC = {"inputFile":["TransmissionCurveFiles/Spitzer_IRAC.I1.dat",
    #                               "TransmissionCurveFiles/Spitzer_IRAC.I2.dat"]}
    # modelIRAC = ObservationModel(zRange,transDictIRAC,spectraDict)
    # modelIRAC.RestSpectrum.addEmissionLine(equivalentWidth=67,
    #                                     lineCenter=4960.3,lineName='OIII-I')
    # modelIRAC.RestSpectrum.addEmissionLine(equivalentWidth=200,
    #                                     lineCenter=5008.2,lineName="OIII-II")
    # modelIRAC.makeObservedSpectra()
    # #### testing 

    # # _ = modelIRAC.RestSpectrum.showSpectrum()
    # # _ = modelIRAC.ObservedSpectra['5.6'].showSpectrum()


    # ####


    # # runAnimation(modelIRAC)

    # modelIRAC.runObservations()
    
    # colourDictIRAC = {}
    # for key,mag in modelIRAC.magnitudes.items():
    #     yminusY = mag[0] - mag[1]
    #     colourDictIRAC[key] = yminusY

    # fig,axs = plt.subplots()
  
    # axs.plot(modelIRAC.zRange,colourDictIRAC.values(),color='r')
    
    # axs.set_xlabel("z")
    # axs.set_ylabel(r"3.6$\mu$m - 4.5$\mu$m Colour")
    # axs.set_xlim(xmin=5,xmax=9)
    # # axs.set_ylim(ymin=-2,ymax=2)
    # # axs.set_ylim(ymin=-0.25,ymax=0.25)

    # axs.axhline(0.0,color='k',ls=":")
    # axs.legend(fontsize=8)

    # plt.show()

if __name__ == '__main__':
    main()
