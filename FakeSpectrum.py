import numpy as np
import matplotlib.pyplot as plt

class FakeSpectrum(object):

    def __init__(self,nuIndex=0, lowerLambda=100,upperLambda=30000,
                 normalisationMag=25,activateIGM=True):
        """
        Instance with attributes wavelength & flux
        """
        self.cSpeed = 3.0E+8
        self.normalisationWavelength = 25000.0 #2.5microns

        self.nuIndex = nuIndex
        self.nuNorm = self.setNuNorm(normalisationMag)

        self.lambdaIndex = nuIndex-2.0
        self.lambdaNorm = self.setLambdaNorm(normalisationMag)

        self.wavelength = np.arange(lowerLambda,upperLambda,10)
        self.frequency = self.fluxdensityUnitSwap(self.wavelength)

        self.f_flux = self.spectrumNuModel(self.frequency,self.nuIndex)
        self.w_flux = self.spectrumLambdaModel(self.wavelength,
                                               self.lambdaIndex)

        if activateIGM:
            self.IGM_absorption()

    def IGM_absorption(self):
        #Flux should be set to 0 at lambda<1.216um or freq>2.4653985e+14Hz
        IGM_bool = (self.wavelength<1216)
        self.f_flux[IGM_bool] = 9.1201084e-60
        self.w_flux[IGM_bool] = (9.1201082e-50 * self.cSpeed
                                /self.wavelength[IGM_bool]**2.0)

    def setLambdaNorm(self,normalisationMag):
        w_flux = self.fluxdensityUnitSwap( self.ABTof_flux(normalisationMag),
                                    'frequency', self.normalisationWavelength )
        return w_flux / self.normalisationWavelength**self.lambdaIndex

    def setNuNorm(self,normalisationMag): #3.6307805E-30
        self.normalisationNu = 1E10*self.cSpeed/self.normalisationWavelength
        return (self.ABTof_flux(normalisationMag)
                / self.normalisationNu**self.nuIndex )

    def fluxdensityUnitSwap(self,fluxdensity,currentUnit='wavelength',
                            wavelength=None):
        """
        Converts from F_lam [erg/s/cm2/A] <-> F_nu [erg/s/cm2/Hz]
        """
        if wavelength == None:
            wavelength = self.wavelength

        if currentUnit == 'wavelength':  #fv=3.34x10^-19 * l**2 * fl
            return fluxdensity * wavelength**2.0 * 1.0E-10 / self.cSpeed
        elif currentUnit == 'frequency': #fl=3x10^18 * fv / l**2
            return fluxdensity * 1.0E+10 * self.cSpeed / wavelength**2.0

    def spectrumNuModel(self,frequency,nuIndex):
        """
        Spectrum modelled as fv=a*frequency**nuIndex
        """
        return self.nuNorm*frequency**nuIndex

    def spectrumLambdaModel(self,wavelength,lambdaIndex):
        """
        Spectrum modelled as fl=b*wavelength**alpha, where alpha=nuIndex-2
        """
        return self.lambdaNorm*wavelength**lambdaIndex

    def f_fluxToAB(self,f_flux):
        return -2.5*np.log10(f_flux)-48.6

    def w_fluxToAB(self,w_flux):
        return -2.5*np.log10(self.fluxdensityUnitSwap(w_flux))-48.6

    def ABTof_flux(self,mag):
        return 10.0**(-(mag+48.6)/2.5)

    def showSpectrum(self):
        _,axs = plt.subplots()
        axs.plot(self.wavelength/1E4,self.f_fluxToAB(self.f_flux))
        axs.set_xlabel(r'$\lambda/\mu$m')
        axs.set_ylabel(r'$m_{AB}$')
        axs.set_ylim(ymax=30,ymin=22)
        axs.invert_yaxis()
        plt.show()

def main():
    nuIndex = 0.0
    sed = FakeSpectrum(nuIndex)
    sed.showSpectrum()

if __name__ == '__main__':
    main()
