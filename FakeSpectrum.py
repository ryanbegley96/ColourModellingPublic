import numpy as np
import matplotlib.pyplot as plt

class FakeSpectrum(object):

    def __init__(self,nuIndex, lowerLambda=100,upperLambda=30000):
        """
        Instance with attributes wavelength & flux
        """
        self.cSpeed = 3.0E+8
        self.nuIndex = nuIndex
        self.lambdaIndex = nuIndex-2.0

        self.wavelength = np.linspace(lowerLambda,upperLambda,100)
        self.frequency = self.fluxdensityUnitSwap(self.wavelength)

        self.f_flux = self.spectrumNuModel(self.frequency,self.nuIndex)
        self.w_flux = self.spectrumLambdaModel(self.wavelength,
                                               self.lambdaIndex)
        self.IGM_absorption()

    def IGM_absorption(self):
        #Flux should be set to 0 at lambda<1.216um or freq>2.4653985e+14Hz
        IGM_bool = (self.wavelength<1216)
        self.f_flux[IGM_bool] = 9.1201084e-60
        self.w_flux[IGM_bool] = (9.1201082e-50 * self.cSpeed
                                /self.wavelength[IGM_bool]**2.0)

    def fluxdensityUnitSwap(self,fluxdensity,currentUnit='wavelength'):
        """
        Converts from F_lam [erg/s/cm2/A] <-> F_nu [erg/s/cm2/Hz]
        """
        if currentUnit == 'wavelength':  #fv=3.34x10^-19 * l**2 * fl
            return fluxdensity * self.wavelength**2.0 * 1.0E-10 / self.cSpeed
        elif currentUnit == 'frequency': #fl=3x10^18 * fv / l**2
            return fluxdensity * 1.0E+10 * self.cSpeed / self.wavelength**2.0

    def spectrumNuModel(self,frequency,nuIndex):
        """
        Spectrum modelled as fv=a*frequency**nuIndex
        """
        self.nuNorm = 3.6307805E-30
        return self.nuNorm*frequency**nuIndex

    def spectrumLambdaModel(self,wavelength,lambdaIndex):
        """
        Spectrum modelled as fl=b*wavelength**alpha, where alpha=nuIndex-2
        """
        self.lambdaNorm = 1.08537E-11
        return self.lambdaNorm*wavelength**lambdaIndex

    def f_fluxToAB(self,f_flux):
        return -2.5*np.log10(f_flux)-48.6

    def w_fluxToAB(self,w_flux):
        return -2.5*np.log10(self.fluxdensityUnitSwap(w_flux))-48.6


    def spectrumMagPlotter(self):
        fig,axs = plt.subplots()
        axs.plot(self.wavelength/1E4,self.f_fluxToAB(self.f_flux))
        axs.set_xlabel(r'$\lambda/\mu$m')
        axs.set_ylabel(r'$m_{AB}$')
        axs.set_ylim(ymax=30,ymin=22)
        axs.invert_yaxis()
        plt.show()

def main():
    nuIndex = 0.0
    sed = FakeSpectrum(nuIndex)
    sed.spectrumMagPlotter()

if __name__ == '__main__':
    main()
