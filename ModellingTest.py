import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

from TransmissionCurve import TransmissionCurve
from FakeSpectrum import FakeSpectrum
from MagnitudeModel import MagnitudeModel
from ObservationModel import ObservationModel

def main():
    lineType = 'Gaussian'
    spectraDict = {"nuIndex":0.0,
                   "lowerLambda":100,
                   "upperLambda":30000,
                   "normalisationMag":25,
                   "activateIGM":True,
                   "activateLya":lineType}

    transDict = {"inputFile":["TransmissionCurveFiles/Subaru_HSC.Y.dat",
                              "TransmissionCurveFiles/Paranal_VISTA.Y.dat"]}

    zRange = [0,10]
    model = ObservationModel(zRange,transDict,spectraDict)

    model.runObservations()

    obj1 = model.ObservedSpectra['0.0']
    filter1 = model.transCurves[0]
    filter2 = model.transCurves[1]

    colourDict = {}
    for key,mag in model.magnitudes.items():
        yminusY = mag[0] - mag[1]
        colourDict[key] = yminusY

    rossCalculations = np.loadtxt("Ross_delta_Y_colour.dat")
    rossZ,rossColour = rossCalculations[:,0],rossCalculations[:,-1]
    diracCalulations = np.loadtxt("dirac_delta_Y_colour.dat")
    gaussianCalulations = np.loadtxt("gaussian_delta_Y_colour.dat")

    ##running small animation showing the redshifted spectrum
    ##moving through the filters
    # plt.ion() #needed if doing manual loop animation
    fig,axs = plt.subplots()
    spectrumPlot, = plt.plot(obj1.wavelength/1E4,obj1.f_fluxToAB(obj1.f_flux),'k')
    axs.set_xlabel(r'$\lambda/\mu$m')
    axs.set_ylabel(r'$m_{AB}$')
    axs.set_xlim(xmin=0.,xmax=3.0)
    axs.set_ylim(ymax=30,ymin=22)
    axs.invert_yaxis()

    props = dict(boxstyle='Square', alpha=0.0)
    txtBox = axs.text(0.025, 0.975, "z="+str(zRange[0]), transform=axs.transAxes, fontsize=12,
    verticalalignment='top', bbox=props)

    ax2 = axs.twinx()
    ax2.fill_between(filter1.wavelength/1E4,filter1.transmission,edgecolor='r',facecolor='r',alpha=0.5,label="HSC-Y")
    ax2.fill_between(filter2.wavelength/1E4,filter2.transmission,edgecolor='b',facecolor='b',alpha=0.5,label="UVISTA-Y")
    ax2.set_ylabel("Transmission/%")
    ax2.set_ylim(ymin=0,ymax=1)
    ax2.legend(fontsize='10')

    
    axins = inset_axes(axs, width=1.9, height=1.4, bbox_to_anchor=(0.99,0.485),
                        bbox_transform=axs.transAxes)
    axins.plot(model.zRange,colourDict.values(),color='k')
    axins.plot(rossZ,rossColour,label="Ross' Calc.",color='r')
    axins.plot(diracCalulations[:,0],diracCalulations[:,1],color='b',label="Dirac")
    axins.plot(gaussianCalulations[:,0],gaussianCalulations[:,1],color='g',label="Gaussian")
    axins.axhline(0.0,color='gray',ls="--")    

    currentColour, = axins.plot(0.0,0.0,'k',marker='o',markersize=5)
    axins.set_xlabel(r"z", fontsize=7.5, labelpad=0.0005)
    axins.set_ylabel(r"y-Y Colour", fontsize=7.5, labelpad=0.1)
    axins.xaxis.set_tick_params(labelsize=7.5)
    axins.yaxis.set_tick_params(labelsize=7.5)
    axins.set_xlim(xmin=6,xmax=7.5)
    axins.set_ylim(ymin=-2,ymax=2)
    axins.legend(fontsize=6)

    # for idx,spectrumTuple in enumerate(model.ObservedSpectra.items()):
    #     plt.pause(0.1)
    #     plot1.set_xdata(spectrumTuple[1].wavelength/1E4)
    #     txtBox.set_text("z="+str(np.round(model.zRange[idx],1)))
    #     fig.canvas.draw()

    def ColourAnim(zKey):
        zKey = np.round(zKey,1)
        #Takes the spectrumTuple of ('z',model.object) & updates plot
        spectrum = model.ObservedSpectra[str(zKey)]
        colourValue = colourDict[str(zKey)]
            
        spectrumPlot.set_xdata(spectrum.wavelength/1E4)
        txtBox.set_text("z="+str(zKey))
        currentColour.set_data(zKey,colourValue)

        return spectrumPlot,txtBox,currentColour

    anim = FuncAnimation(fig,ColourAnim,frames=model.zRange[56:84],blit=True)
    plt.show()
    # anim.save('test.mp4')

    #For the y-Y colour vs z plot
    # fig,axs = plt.subplots()
    # axs.plot(model.zRange,colourDict.values(),label=lineType)
    # axs.set_xlabel("z")
    # axs.set_ylabel("y-Y Colour")
    # axs.set_xlim(xmin=6,xmax=7.5)
    # axs.set_ylim(ymin=-2,ymax=2)

    # axs.plot(rossZ,rossColour,label="Ross' Calc.")
    # axs.plot(diracCalulations[:,0],diracCalulations[:,1],label="Dirac")
    # axs.axhline(0.0,color='k',ls="--")
    # axs.legend(fontsize=8)

    # plt.show()

if __name__ == '__main__':
    main()
