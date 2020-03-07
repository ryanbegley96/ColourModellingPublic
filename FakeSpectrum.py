import numpy as np
import matplotlib.pyplot as plt

class FakeSpectrum(object):

    def __init__(self,nuIndex, lowerLambda=0.1,upperLambda=3.0):
        """
        Instance with attributes wavelength & flux
        """
        self.normalisation = 3.6307805 #e-30
        self.wavelength = np.linspace(lowerLambda,upperLambda,100)
        self.flux = self.spectrumNuModel(self.wavelength,nuIndex)

    def spectrumNuModel(self,wavelength,nu):
        """
        Spectrum modelled as f=a*lambda**nu w/ full IGM absorption.
        Normalised to ABmag=25 w/ ZP=48.6 -> cgs units [erg/s^1/cm^2/Hz]
        a = 3.6307805e-30 #removing e-30 for now
        """
        flux = self.normalisation*self.wavelength**nu
        flux[self.wavelength < 1.215] = 0.0
        
        return flux

    def spectrumPlotter(self):
        fig,axs = plt.subplots()
        axs.plot(self.wavelength,self.flux)
        plt.show()

def main():
    nuIndex = 0.0
    sed = FakeSpectrum(nuIndex)
    sed.spectrumPlotter()

if __name__ == '__main__':
    main()