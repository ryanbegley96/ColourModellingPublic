import numpy as np
import matplotlib.pyplot as plt

from TransmissionCurve import TransmissionCurve
from FakeSpectrum import FakeSpectrum

class MagnitudeModel(object):

    def __init__(self,spectaDict,transmissionDict):

        self.spectrum = FakeSpectrum(**spectaDict)
        self.spectrumAB = FakeSpectrum(0.0,1E2,3E4,0.0,False)

        self.transCurve = TransmissionCurve(**transmissionDict)
        self.transFunction = self.transCurve.returnInterpolationFunc()
        self.resampleTransmissionCurve()

    def resampleTransmissionCurve(self):
        """
        Sample Tranmission Curve for values of wavelength within the defined
        range and zero otherwise.
        """
        samplingRange = (np.min(self.transCurve.wavelength),
                         np.max(self.transCurve.wavelength))
        samplingWavelengthBool = ((self.spectrum.wavelength>samplingRange[0])
                                  & (self.spectrum.wavelength<samplingRange[1]))

        #produce arrays of value over the range of the filter.
        self.samplingWavelength = self.spectrum.wavelength[samplingWavelengthBool]
        self.sampledTransCurve = self.transFunction(self.samplingWavelength)
        self.sampledFlux = self.spectrum.w_flux[samplingWavelengthBool]
        self.sampledFluxAB = self.spectrumAB.w_flux[samplingWavelengthBool]

    def computeMagnitude(self):
        print(np.sum(self.sampledFlux*self.samplingWavelength*self.sampledTransCurve)/
              np.sum(self.sampledFluxAB*self.samplingWavelength*self.sampledTransCurve))


def main():
    spectraDict = {"nuIndex":0.0,
                   "lowerLambda":100,
                   "upperLambda":30000,
                   "normalisationMag":25,
                   "activateIGM":True}

    transDict = {"inputFile":"TransmissionCurveFiles/Subaru_HSC.g_filter.dat"}
    model = MagnitudeModel(spectraDict,transDict)
    model.computeMagnitude()

if __name__ == '__main__':
    main()
