import copy
import time
import numpy as np
import matplotlib.pyplot as plt

from TransmissionCurve import TransmissionCurve
from FakeSpectrum import FakeSpectrum
from MagnitudeModel import MagnitudeModel

class ObservationModel():

    def __init__(self,zRange,transDict,spectraDict):

        self.zRange = np.arange(*zRange,0.1)
        self.transCurves = [TransmissionCurve(file)
                            for file in transDict['inputFile']]
        self.spectrumAB = FakeSpectrum(0.0,1E2,3E4,0.0,False)
        self.RestSpectrum = FakeSpectrum(**spectraDict)
        self.ObservedSpectra = self.makeObservedSpectra()
        
    def makeObservedSpectra(self):
        ObservedSpectra = {}
        for z in self.zRange:
            observedSpectrum = copy.deepcopy(self.RestSpectrum)
            observedSpectrum.wavelength = observedSpectrum.wavelength * (1.0+z)
            ObservedSpectra[str(z)] = observedSpectrum
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
    print("mag :", model.magnitudes)

    obj1 = model.ObservedSpectra['0.0']
    filter1 = model.transCurves[0]
    filter2 = model.transCurves[1]

    # plt.ion()
    # fig,axs = plt.subplots()
    # plot1, = plt.plot(obj1.wavelength/1E4,obj1.f_fluxToAB(obj1.f_flux),'k')
    # axs.set_xlabel(r'$\lambda/\mu$m')
    # axs.set_ylabel(r'$m_{AB}$')
    # axs.set_xlim(xmin=0.,xmax=3.0)
    # axs.set_ylim(ymax=30,ymin=22)
    # axs.invert_yaxis()

    # ax2 = axs.twinx()
    # ax2.fill_between(filter1.wavelength/1E4,filter1.transmission,edgecolor='r',facecolor='r',alpha=0.5,label="HSC-Y")
    # ax2.fill_between(filter2.wavelength/1E4,filter2.transmission,edgecolor='b',facecolor='b',alpha=0.5,label="UVISTA-Y")
    # ax2.set_ylabel("Transmission/%")
    # ax2.set_ylim(ymin=0,ymax=1)
    # ax2.legend(fontsize='10')

    # for _,spectrum in model.ObservedSpectra.items():
    #     plt.pause(0.1)
    #     plot1.set_xdata(spectrum.wavelength/1E4)
    #     fig.canvas.draw()

    # colourList = 
    # fig,axs = plt.subplots()


    plt.show()

if __name__ == '__main__':
    main()

"""
def test(spectraDict,transDict):
    model = ObservationModel([0,1],{"inputFile":["testTransmissionCurve.txt"]},spectraDict)
    magModel = MagnitudeModel(model.RestSpectrum, model.spectrumAB, 
                              model.transCurves[0])
    print("calc : ",magModel.magnitude)
    # print("sampled vals : ",magModel.sampledFlux)
    
    fig,axs = plt.subplots()
    plt.plot(model.RestSpectrum.wavelength/1E4,model.RestSpectrum.f_fluxToAB(model.RestSpectrum.f_flux))
    plt.plot(model.spectrumAB.wavelength/1E4,model.spectrumAB.f_fluxToAB(model.spectrumAB.f_flux))
    axs.set_xlabel(r'$\lambda/\mu$m')
    axs.set_ylabel(r'$m_{AB}$')
    axs.set_xlim(xmin=0.,xmax=3.0)
    # axs.set_ylim(ymax=30,ymin=22)
    filter1 = model.transCurves[0]
    axs.invert_yaxis()
    ax2 = axs.twinx()
    # ax2.plot(filter1.wavelength/1E4,filter1.transmission,c='r')
    ax2.fill_between(filter1.wavelength/1E4,filter1.transmission,edgecolor='r',facecolor='r',alpha=0.5)
    ax2.set_ylabel("Transmission/%")
    ax2.set_ylim(ymin=0,ymax=1)
"""