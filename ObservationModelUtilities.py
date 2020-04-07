"""
File containing useful example functions for use with ObservationModel
and other classes. Particulary called into ModellingTest.py.
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

from ObservationModel import ObservationModel

def exampleModel(zRange,transDict,spectraDict,*comparisonsToShow):
    """Exampel usage & showing results.
    Keyword arguments:
    zRange -- (min,max) to run model over
    transDict -- dict with key "inputFile" containing list of filenames.
    spectraDict -- dict with keys as keywords setting Spectrum Class.
    *comparisonsToShow -- args in the form of tuple(label,data,colour) with data
                          being two columns : z,y-Y_colour
    """
    spectraDict["activateLya"]="Gaussian"
    
    model = ObservationModel(zRange,transDict,spectraDict)
    model.runObservations()

    colourDict = {}
    for key,mag in model.magnitudes.items():
        yminusY = mag[0] - mag[1]
        colourDict[key] = yminusY

    ###For the y-Y colour vs z plot
    fig,axs = plt.subplots()
  
    axs.plot(model.zRange,colourDict.values(),label='Gaussian, Intr. EW=200',color='r')
    
    axs.set_xlabel("z")
    axs.set_ylabel("y-Y Colour")
    axs.set_xlim(xmin=6,xmax=7.5)
    axs.set_ylim(ymin=-2,ymax=2)

    for comparison in comparisonsToShow:
        label,data,colour = comparison
        axs.plot(data[:,0],data[:,-1],label="Comparison:"+label,c=colour)
    
    axs.axhline(0.0,color='k',ls=":")
    axs.legend(fontsize=8)

    plt.show()


def objectOverlayPlotter(axs, objectData, plottingType="Default"):
    """ Plots the object redshift vs colour from object data onto axs.
    Keyword arguments:
    axs -- matplotlib figure axs to plot object data on.
    objectData -- array containing object data with one object per row
                  and columns: ID,z,...,colour.
    plottingType -- options for how to plot objects:Default,All,None.
    """
    ##plotting for the objects
    if plottingType is "Default":
        print("Plotting default manually selected Objects.")

        minColourIdx = np.argmin(objectData[:,-1])
        maxColourIdx = np.argmax(objectData[:,-1])

        objZ = objectData[minColourIdx,1]
        objC = objectData[minColourIdx,-1]
        axs.plot(objZ,objC,markersize=8,marker='*',color='k')
        props = dict(boxstyle='Square', alpha=0.0)
        txtBox = axs.text(objZ+0.02,objC+0.05, "ID:"+str(int(objectData[minColourIdx,0])),fontsize=8, bbox=props)

        objZ1454270 = objectData[3,1]
        objC1454270 = objectData[3,-1]
        axs.plot(objZ1454270,objC1454270,markersize=8,marker='*',color='k')
        props = dict(boxstyle='Square', alpha=0.0)
        txtBox = axs.text(objZ1454270+0.02,objC1454270+0.05, "ID:"+str(int(objectData[2,0])),fontsize=8, bbox=props)
        
        objZ22022334 = objectData[8,1]
        objC22022334= objectData[8,-1]
        axs.plot(objZ22022334,objC22022334,markersize=8,marker='*',color='k')
        props = dict(boxstyle='Square', alpha=0.0)
        txtBox = axs.text(objZ22022334+0.02,objC22022334+0.05, "ID:"+str(int(objectData[8,0])),fontsize=8, bbox=props)

    elif plottingType is "All":
        print("Plotting all Objects.")

        props = dict(boxstyle='Square', alpha=0.0)
        
        for obj in objectData:
            objZ = obj[1]
            objC = obj[-1]
            axs.plot(objZ,objC,markersize=8,marker='*',color='k')
            txtBox = axs.text(objZ+0.02,objC+0.05, "ID:"+str(int(obj[0])),
                            fontsize=8, bbox=props)
    else:
        print("No objects overlayed onto model plot.")


def setupAnimationFigure(model,colourDict):
    """ Setup figure required for animation including
        artists that will be edited.
    Keyword arguments:
    model -- model initialisation including spectra & filter to run over z.
    colourDict -- results of running the modelObservations, dict of form
                  {key=str(redshift,1dp):value=colourvalue}
    """
    restObj = model.RestSpectrum
    filter1 = model.transCurves[0]
    filter2 = model.transCurves[1]
    
    rossCalculations = np.loadtxt("Ross_delta_Y_colour.dat")
    rossCalculations_v2 = np.loadtxt("Ross_delta_Y_colour_v2.dat")
    diracCalulations = np.loadtxt("dirac_delta_Y_colour.dat")
    gaussianCalulations = np.loadtxt("gaussian_delta_Y_colour.dat")
    
    fig,axs = plt.subplots()
    spectrumPlot, = plt.plot(restObj.wavelength/1E4,restObj.f_fluxToAB(restObj.f_flux),'k')
    axs.set_xlabel(r'$\lambda/\mu$m')
    axs.set_ylabel(r'$m_{AB}$')
    axs.set_xlim(xmin=0.,xmax=3.0)
    axs.set_ylim(ymax=30,ymin=22)
    axs.invert_yaxis()

    props = dict(boxstyle='Square', alpha=0.0)
    txtBox = axs.text(0.025, 0.975, "z="+str(model.zRange[0]), 
            transform=axs.transAxes, fontsize=12,verticalalignment='top', 
            bbox=props)

    ax2 = axs.twinx()
    ax2.fill_between(filter1.wavelength/1E4,filter1.transmission,edgecolor='r',facecolor='r',alpha=0.5,label="HSC-Y")
    ax2.fill_between(filter2.wavelength/1E4,filter2.transmission,edgecolor='b',facecolor='b',alpha=0.5,label="UVISTA-Y")
    ax2.set_ylabel("Transmission/%")
    ax2.set_ylim(ymin=0,ymax=1)
    ax2.legend(fontsize='10')

    
    axins = inset_axes(axs, width=1.9, height=1.4, bbox_to_anchor=(0.99,0.485),
                        bbox_transform=axs.transAxes)
    axins.plot(model.zRange,colourDict.values(),color='k')
    axins.plot(rossCalculations_v2[:,0],rossCalculations_v2[:,-1],label="Ross' Calc. v2",color='r')
    axins.plot(diracCalulations[:,0],diracCalulations[:,-1],color='b',label="Dirac")
    axins.plot(gaussianCalulations[:,0],gaussianCalulations[:,-1],color='g',label="Gaussian")
    axins.axhline(0.0,color='gray',ls="--")    

    currentColour, = axins.plot(0.0,0.0,'k',marker='o',markersize=5)
    axins.set_xlabel(r"z", fontsize=7.5, labelpad=0.0005)
    axins.set_ylabel(r"y-Y Colour", fontsize=7.5, labelpad=0.1)
    axins.xaxis.set_tick_params(labelsize=7.5)
    axins.yaxis.set_tick_params(labelsize=7.5)
    axins.set_xlim(xmin=6,xmax=7.5)
    axins.set_ylim(ymin=-2,ymax=2)
    axins.legend(fontsize=6)
    return fig,axs,spectrumPlot,txtBox,currentColour

def runAnimation(model):
    """ Run animation demonstrating how model obseration performed.
    Keyword arguments:
    model --model initialisation including spectra & filter to run over z.
    """
    model.runObservations()

    colourDict = {}
    for key,mag in model.magnitudes.items():
        yminusY = mag[0] - mag[1]
        colourDict[key] = yminusY

    fig,axs,spectrumPlot,txtBox,currentColour = setupAnimationFigure(model,colourDict)

    def ModelAnimUpdater(zKey,model,colourDict):
        zKey = np.round(zKey,1)
        #Takes the spectrumTuple of ('z',model.object) & updates plot
        spectrum = model.ObservedSpectra[str(zKey)]
        colourValue = colourDict[str(zKey)]
            
        spectrumPlot.set_xdata(spectrum.wavelength/1E4)
        txtBox.set_text("z="+str(zKey))
        currentColour.set_data(zKey,colourValue)

        return spectrumPlot,txtBox,currentColour

    anim = FuncAnimation(fig,ModelAnimUpdater,frames=model.zRange[56:84],
                        fargs=(model,colourDict),blit=True)
    plt.show()