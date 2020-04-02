import copy

import numpy as np
# import matplotlib.pyplot as plt
# from matplotlib.animation import FuncAnimation
# from mpl_toolkits.axes_grid1.inset_locator import inset_axes

from TransmissionCurve import TransmissionCurve
from MagnitudeModel import MagnitudeModel
# from FakeSpectrum import FakeSpectrum
from Spectrum import BaseSpectrum,FakeSpectrum


class ObservationModel():

    def __init__(self,zRange,transDict,spectraDict):

        self.zRange = np.arange(*zRange,0.1)
        self.transCurves = [TransmissionCurve(file)
                            for file in transDict['inputFile']]
        self.spectrumAB = FakeSpectrum(0.0,1E2,3E4,0.0,False,False)
        self.RestSpectrum = FakeSpectrum(**spectraDict)
        self.makeObservedSpectra()
        # self.ObservedSpectra = self.makeObservedSpectra()

    def makeObservedSpectra(self):
        ObservedSpectra = {}
        for z in self.zRange:
            observedSpectrum = copy.deepcopy(self.RestSpectrum)
            observedSpectrum.wavelength = observedSpectrum.wavelength * (1.0+z)
            ObservedSpectra[str(np.round(z,1))] = observedSpectrum
        self.ObservedSpectra = ObservedSpectra
        #return ObservedSpectra

    def runObservations(self):
        self.magnitudes = {}
        for zKey,spectrum in self.ObservedSpectra.items():
            magList = []
            for transCurve in self.transCurves:
                magModel = MagnitudeModel(spectrum,self.spectrumAB,transCurve)
                magList.append(magModel.magnitude)
            self.magnitudes[zKey] = magList
