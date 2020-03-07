import numpy as np
import matplotlib.pyplot as plt

class FakeSpectrum(object):

    def __init__(self,nuIndex, lowerLambda=0.1,upperLambda=3.0):
        """
        Instance with attributes wavelength & flux
        """
        # self.normalisation = 3.6307805 #e-30
        self.cSpeed = 3.0E+8
        self.nuIndex = nuIndex
        self.lambdaIndex = nuIndex-2.0

        self.wavelength = np.linspace(lowerLambda,upperLambda,100)
        self.frequency = self.fluxdensityUnitSwap(self.wavelength)#, 'wavelength')

        self.f_flux = self.spectrumNuModel(self.frequency,self.nuIndex)
        self.w_flux = self.spectrumLambdaModel(self.wavelength,self.lambdaIndex)

    def fluxdensityUnitSwap(self,fluxdensity,currentUnit='wavelength'):
        """
        Converts from F_lam [erg/s/cm2/A] <-> F_nu [erg/s/cm2/Hz]
        """
        if currentUnit == 'wavelength':  #fv=3x10^18 * fl / l**2
            return fluxdensity * 1.0E+10 * self.cSpeed / self.wavelength**2.0
        elif currentUnit == 'frequency': #fl=1x10^10 * fv * l**2 / c
            return fluxdensity * self.wavelength**2.0 * 1.0E-10 / self.cSpeed


    def spectrumNuModel(self,frequency,nuIndex):
        """
        Spectrum modelled as fv=a*frequency**nuIndex
        """
        return frequency**nuIndex
        
    def spectrumLambdaModel(self,wavelength,lambdaIndex):
        """
        Spectrum modelled as fl=b*wavelength**alpha, where alpha=nuIndex-2
        """
        return wavelength**lambdaIndex    

    def spectrumPlotter(self):
        fig,axs = plt.subplots()
        axs.plot(self.wavelength,self.flux)
        plt.show()

def main():
    nuIndex = 0.0
    sed = FakeSpectrum(nuIndex)
    # sed.spectrumPlotter()

if __name__ == '__main__':
    main()
