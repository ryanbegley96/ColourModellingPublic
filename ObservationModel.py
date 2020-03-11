import copy
import numpy as np
import matplotlib.pyplot as plt

from TransmissionCurve import TransmissionCurve
from FakeSpectrum import FakeSpectrum
from MagnitudeModel import MagnitudeModel

class ObservationModel():

    def __init__(self,zRange,transDict,spectraDict):

        self.zRange = np.arange(*zRange,4)
        self.transCurves = [TransmissionCurve(file)
                            for file in transDict['inputFile']]
        self.spectrumAB = FakeSpectrum(0.0,1E2,3E4,0.0,False)
        self.RestSpectrum = FakeSpectrum(**spectraDict)
        self.ObservedSpectra = self.makeObservedSpectra()

    def makeObservedSpectra(self):
        ObservedSpectra = []
        for z in self.zRange:
            observedSpectrum = copy.deepcopy(self.RestSpectrum)
            observedSpectrum.wavelength = observedSpectrum.wavelength * (1.0+z)
            ObservedSpectra.append(observedSpectrum)
        return ObservedSpectra

    def runObservations(self):
        self.magnitudes = []
        for spectrum in self.ObservedSpectra:
            magModel = MagnitudeModel(spectrum,self.spectrumAB,self.transCurves[0])
            self.magnitudes.append(magModel.magnitude)

def main():
    spectraDict = {"nuIndex":0.0,
                   "lowerLambda":100,
                   "upperLambda":30000,
                   "normalisationMag":25,
                   "activateIGM":True}

    transDict = {"inputFile":["TransmissionCurveFiles/Subaru_HSC.g_filter.dat",
                              "TransmissionCurveFiles/Paranal_VISTA.Y.dat"]}

    zRange = [0,5]
    model = ObservationModel(zRange,transDict,spectraDict)
    model.runObservations()
    print(model.magnitudes)

    obj1 = model.ObservedSpectra[0]
    obj2 = model.ObservedSpectra[1]
    filter1 = model.transCurves[0]
    filter2 = model.transCurves[1]

    fig,axs = plt.subplots()
    plt.plot(obj1.wavelength/1E4,obj1.f_fluxToAB(obj1.f_flux))
    plt.plot(obj2.wavelength/1E4,obj2.f_fluxToAB(obj2.f_flux))
    axs.set_xlabel(r'$\lambda/\mu$m')
    axs.set_ylabel(r'$m_{AB}$')
    axs.set_xlim(xmin=0.,xmax=3.0)
    axs.set_ylim(ymax=30,ymin=22)
    axs.invert_yaxis()

    ax2 = axs.twinx()
    # ax2.plot(filter1.wavelength/1E4,filter1.transmission,c='r')
    ax2.fill_between(filter1.wavelength/1E4,filter1.transmission,edgecolor='r',facecolor='r',alpha=0.5)
    # ax2.fill_between(filter1.wavelength/1E4,filter1.transmission,edgecolor='r',facecolor='r',alpha=0.5)
    ax2.set_ylabel("Transmission/%")
    ax2.set_ylim(ymin=0,ymax=1)

    plt.show()

if __name__ == '__main__':
    main()
