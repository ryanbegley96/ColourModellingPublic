import numpy as np
import matplotlib.pyplot as plt

class MagnitudeModel(object):

    def __init__(self,spectrum,spectrumAB,transCurve):

        self.spectrum = spectrum
        self.spectrumAB = spectrumAB
        self.transCurve = transCurve

        resampledResults = self.resampleSpectrum()
        self.samplingWavelength = resampledResults["Wavelength"]
        self.sampledFlux = resampledResults["Flux"]
        self.sampledFluxAB = resampledResults["FluxAB"]
        self.sampledTransCurve = self.transCurve.returnInterpolationFunc()(
                                    self.samplingWavelength)

        self.magnitude = self.calculateMagnitude()


    def resampleSpectrum(self):

        samplingRange = (np.min(self.transCurve.wavelength),
                          np.max(self.transCurve.wavelength))
        samplingWavelengthBool = ((self.spectrum.wavelength>samplingRange[0])
                                  & (self.spectrum.wavelength<samplingRange[1]))
        # produce arrays of value over the range of the filter.
        self.samplingWavelength = self.spectrum.wavelength[samplingWavelengthBool]
        self.sampledTransCurve = self.transCurve.returnInterpolationFunc()(
                                                 self.samplingWavelength)
        self.sampledFlux = self.spectrum.w_flux[samplingWavelengthBool]
        self.sampledFluxAB = self.spectrumAB.w_flux[samplingWavelengthBool]
        result = {"Wavelength":self.spectrum.wavelength[samplingWavelengthBool],
                "Flux":self.spectrum.w_flux[samplingWavelengthBool],
                "FluxAB":self.spectrumAB.w_flux[samplingWavelengthBool],
                }
        return result


    def calculateMagnitude(self):
        integ = np.sum(self.sampledFlux*self.samplingWavelength*
            self.sampledTransCurve)
        norm = np.sum(self.sampledFluxAB*self.samplingWavelength*
            self.sampledTransCurve)
        return -2.5*np.log10(integ/norm)
