import copy
import time
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

from TransmissionCurve import TransmissionCurve
from FakeSpectrum import FakeSpectrum
from MagnitudeModel import MagnitudeModel

class ObservationModel():

    def __init__(self,zRange,transDict,spectraDict):

        self.zRange = np.arange(*zRange,0.1)
        self.transCurves = [TransmissionCurve(file)
                            for file in transDict['inputFile']]
        self.spectrumAB = FakeSpectrum(0.0,1E2,3E4,0.0,False,False)
        self.RestSpectrum = FakeSpectrum(**spectraDict)
        self.ObservedSpectra = self.makeObservedSpectra()

    def makeObservedSpectra(self):
        ObservedSpectra = {}
        for z in self.zRange:
            observedSpectrum = copy.deepcopy(self.RestSpectrum)
            observedSpectrum.wavelength = observedSpectrum.wavelength * (1.0+z)
            ObservedSpectra[str(np.round(z,1))] = observedSpectrum
        return ObservedSpectra

    def runObservations(self):
        self.magnitudes = {}
        for zKey,spectrum in self.ObservedSpectra.items():
            magList = []
            for transCurve in self.transCurves:
                magModel = MagnitudeModel(spectrum,self.spectrumAB,transCurve)
                magList.append(magModel.magnitude)
            self.magnitudes[zKey] = magList

def main():
    spectraDict = {"nuIndex":0.0,
                   "lowerLambda":100,
                   "upperLambda":30000,
                   "normalisationMag":25,
                   "activateIGM":True}

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

    ##running small animation showing the redshifted spectrum
    ##moving through the filters
    # plt.ion() #needed if doing manual loop animation
    # fig,axs = plt.subplots()
    # spectrumPlot, = plt.plot(obj1.wavelength/1E4,obj1.f_fluxToAB(obj1.f_flux),'k')
    # axs.set_xlabel(r'$\lambda/\mu$m')
    # axs.set_ylabel(r'$m_{AB}$')
    # axs.set_xlim(xmin=0.,xmax=3.0)
    # axs.set_ylim(ymax=30,ymin=22)
    # axs.invert_yaxis()

    # props = dict(boxstyle='Square', alpha=0.0)
    # txtBox = axs.text(0.025, 0.975, "z="+str(zRange[0]), transform=axs.transAxes, fontsize=12,
    # verticalalignment='top', bbox=props)

    # ax2 = axs.twinx()
    # ax2.fill_between(filter1.wavelength/1E4,filter1.transmission,edgecolor='r',facecolor='r',alpha=0.5,label="HSC-Y")
    # ax2.fill_between(filter2.wavelength/1E4,filter2.transmission,edgecolor='b',facecolor='b',alpha=0.5,label="UVISTA-Y")
    # ax2.set_ylabel("Transmission/%")
    # ax2.set_ylim(ymin=0,ymax=1)
    # ax2.legend(fontsize='10')

    
    # axins = inset_axes(axs, width=1.9, height=1.4, bbox_to_anchor=(0.99,0.485),
    #                     bbox_transform=axs.transAxes)
    # axins.plot(model.zRange,colourDict.values(),'k')
    # currentColour, = axins.plot(0.0,0.0,'r',marker='o',markersize=8)
    # axins.set_xlabel(r"z", fontsize=7.5, labelpad=0.0005)
    # axins.set_ylabel(r"y-Y Colour", fontsize=7.5, labelpad=0.1)
    # axins.xaxis.set_tick_params(labelsize=7.5)
    # axins.yaxis.set_tick_params(labelsize=7.5)
    # axins.set_xlim(xmin=6,xmax=7.5)
    # axins.set_ylim(ymin=-2,ymax=2)

    # # for idx,spectrumTuple in enumerate(model.ObservedSpectra.items()):
    # #     plt.pause(0.1)
    # #     plot1.set_xdata(spectrumTuple[1].wavelength/1E4)
    # #     txtBox.set_text("z="+str(np.round(model.zRange[idx],1)))
    # #     fig.canvas.draw()

    # def ColourAnim(zKey):
    #     zKey = np.round(zKey,1)
    #     #Takes the spectrumTuple of ('z',model.object) & updates plot
    #     spectrum = model.ObservedSpectra[str(zKey)]
    #     colourValue = colourDict[str(zKey)]
            
    #     spectrumPlot.set_xdata(spectrum.wavelength/1E4)
    #     txtBox.set_text("z="+str(zKey))
    #     currentColour.set_data(zKey,colourValue)

    #     return spectrumPlot,txtBox,currentColour

    # anim = FuncAnimation(fig,ColourAnim,frames=model.zRange[56:84],blit=True)
    # plt.show()
    # anim.save('test.mp4')

    ##For the y-Y colour vs z plot
    # fig,axs = plt.subplots()
    # axs.plot(model.zRange,colourDict.values())
    # axs.set_xlabel("z")
    # axs.set_ylabel("y-Y Colour")
    # axs.set_xlim(xmin=6,xmax=7.5)
    # axs.set_ylim(ymin=-2,ymax=2)

    # plt.show()

if __name__ == '__main__':
    main()
