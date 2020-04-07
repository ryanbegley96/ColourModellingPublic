import copy 
import numpy as np
import matplotlib.pyplot as plt

class BaseSpectrum(object):

    def __init__(self,nuIndex,lowerLambda,upperLambda,normalisationMag):
        
        self.cSpeed = 3.0E+8
        self.normalisationMag = normalisationMag
        self.normalisationWavelength = 25000.0 #2.5microns
        self.lambdaSpacing = 1

        self.nuIndex = nuIndex
        self.nuNorm = self.setNuNorm(normalisationMag)

        self.lambdaIndex = nuIndex-2.0
        self.lambdaNorm = self.setLambdaNorm(normalisationMag)

        self.wavelength = np.arange(lowerLambda,upperLambda,self.lambdaSpacing)
        self.frequency = 1.0E-10 * self.cSpeed /self.wavelength

        self.f_flux = self.spectrumNuModel(self.frequency,self.nuIndex)
        self.w_flux = self.spectrumLambdaModel(self.wavelength,
                                               self.lambdaIndex)

    def setNuNorm(self,normalisationMag): #3.6307805E-30
        self.normalisationNu = 1E10*self.cSpeed/self.normalisationWavelength
        return (self.MagToFlux(normalisationMag)
                / self.normalisationNu**self.nuIndex )

    def setLambdaNorm(self,normalisationMag):
        w_flux = self.fluxdensityUnitSwap( self.MagToFlux(normalisationMag),
                                    'frequency', self.normalisationWavelength )
        return w_flux / self.normalisationWavelength**self.lambdaIndex

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

    def fluxdensityUnitSwap(self,fluxdensity,currentUnit='wavelength',
                        wavelength=None):
        """
        Converts from F_lam [erg/s/cm2/A] <-> F_nu [erg/s/cm2/Hz]
        """
        if wavelength is None:
            wavelength = self.wavelength

        if currentUnit == 'wavelength':  #fv=3.34x10^-19 * l**2 * fl
            return fluxdensity * wavelength**2.0 * 1.0E-10 / self.cSpeed
        elif currentUnit == 'frequency': #fl=3x10^18 * fv / l**2
            return fluxdensity * 1.0E+10 * self.cSpeed / wavelength**2.0
        else:
            raise ValueError( ("Bad Argument in fluxdensityUnitSwap, currentUnit must"
                            "be 'wavelength' or 'frequency'.")
                        )    

    def f_fluxToAB(self,f_flux):
        return -2.5*np.log10(f_flux)-48.6
    
    def w_fluxToAB(self,w_flux,wavelength=None):
        if wavelength is None:
            return -2.5*np.log10(self.fluxdensityUnitSwap(w_flux))-48.6
        else:
            return -2.5*np.log10(self.fluxdensityUnitSwap(w_flux,'wavelength',
                                wavelength) )-48.6

    def MagToFlux(self,mag):
        return 10.0**(-(mag+48.6)/2.5)

    def showSpectrum(self):
        _,axs = plt.subplots()
        axs.plot(self.wavelength/1E4,self.f_fluxToAB(self.f_flux))
        axs.set_xlabel(r'$\lambda/\mu$m')
        axs.set_ylabel(r'$m_{AB}$')
        axs.set_xlim(xmin=0.,xmax=3.0)
        axs.set_ylim(ymax=30,ymin=22)
        axs.invert_yaxis()
        plt.show()
        return axs

class FakeSpectrum(BaseSpectrum):

    def __init__(self,nuIndex=0, lowerLambda=100,upperLambda=30000,
                 normalisationMag=25,activateIGM=True,activateLya=None):
        super().__init__(nuIndex,lowerLambda,upperLambda,normalisationMag)
        
        self.emissionLineRecord = {}
        self.emissionLineData = {}
        self.LyaType = activateLya
        self.IGM_absorptionBool = activateIGM 

        if activateLya is "Dirac":
            self.implementDiracLya()
              
        if activateLya is "Gaussian":
            self.implementGaussianLya()

        if activateIGM:
            self.IGM_absorption()

    
    def addEmissionLine(self,equivalentWidth=200,lineCenter=1215.7,
                        lineName='lya'):
        """ Add emission line to the fake spectrum.
        
        Keyword arguments:
        equivalentWidth -- intrinsic EW of line
        lineCenter -- wavelength of line
        lineName -- string name of line
        """
        lineWidth = 5.0 #default value
        integLimits = (lineCenter-5.*lineWidth,lineCenter+5.*lineWidth)
        
        idx = np.argmin(np.abs(self.wavelength-lineCenter))
        lineConst = equivalentWidth * self.w_flux[idx] 

        wavelengthBool = ( (self.wavelength > integLimits[0]) & 
                            (self.wavelength < integLimits[1]) )
        
        lineWavelengths = self.wavelength[wavelengthBool]
        continuumFlux = copy.deepcopy(self.w_flux)[wavelengthBool]

        lineFlux = self.gaussianFunc(lineCenter,lineWidth,lineConst,
                                     lineWavelengths)
        
        if self.IGM_absorptionBool:
                lineFlux[(lineWavelengths<1216)] = 0.0

        self.w_flux[wavelengthBool] += lineFlux
        self.f_flux[wavelengthBool] += self.fluxdensityUnitSwap(lineFlux,
                                            'wavelength',lineWavelengths)
        equivalentWidthCalc = self.equivWidthIntegral(lineFlux,continuumFlux,
                                                    self.lambdaSpacing)
    
        self.emissionLineRecord[lineName] = (lineCenter,equivalentWidthCalc)
        self.emissionLineData[lineName] = (lineWavelengths,lineFlux)

    def removeEmissionLine(self,lineName):
        """ Remove emission line previously added by addEmissionLine().
        Uses the emissionLineData attribute dict storing line wl,f.
        lineName must be made in emissionLineRecord
        Keyword Arguments:
        lineName -- string name of line
        """
        if lineName in self.emissionLineRecord:

            lineWavelengths,lineFlux = self.emissionLineData[lineName]
            wavelengthBool = np.isin(self.wavelength,lineWavelengths)
            
            self.w_flux[wavelengthBool]-=lineFlux
            
            self.f_flux[wavelengthBool]-=self.fluxdensityUnitSwap(lineFlux,
                                            'wavelength',lineWavelengths)
            
            del self.emissionLineRecord[lineName]
            del self.emissionLineData[lineName]
        
        else:
            raise ValueError("No Line present in the fake spectrum called:"
                            +lineName)    
            #If no code break to be implemented call a pass instead.

    def implementGaussianLya(self):
        equivalentWidth = 400.0
        lineCenter = 1215.7 #wavelength of Lya emission line
        lineWidth = 5.0 #sigma of Gaussian in Angstroms
        integLimits = (lineCenter-5.*lineWidth,lineCenter+5.*lineWidth)
        
        idx = np.argmin(np.abs(self.wavelength-1215.7)) #fc val @ 1216A
        lineConst = equivalentWidth * self.w_flux[idx] #7.366E-18#erg/s/cm2/A

        wavelengthBool = ( (self.wavelength > integLimits[0]) & 
                            (self.wavelength < integLimits[1]) )

        lineWavelengths = self.wavelength[wavelengthBool]
        continuumFlux = copy.deepcopy(self.w_flux)[wavelengthBool]

        lineFlux = self.gaussianFunc(lineCenter,lineWidth,lineConst,
                                     lineWavelengths)
        self.w_flux[wavelengthBool] += lineFlux
        self.f_flux[wavelengthBool] += self.fluxdensityUnitSwap(lineFlux,
                                            'wavelength',lineWavelengths)
        equivalentWidthCalc = self.equivWidthIntegral(lineFlux,continuumFlux,
                                                    self.lambdaSpacing)
        self.emissionLineRecord['lya'] = (lineCenter,equivalentWidthCalc)

    def implementDiracLya(self):
        equivalentWidth = 200.0
        idx = np.argmin(np.abs(self.wavelength-1215.7))
        A = equivalentWidth * self.w_flux[idx] #7.366E-18

        self.w_flux[idx] += A
        self.f_flux[idx] += self.fluxdensityUnitSwap(A,'wavelength',1216.0)

    def IGM_absorption(self):
        #Flux should be set to 0 at lambda<1.216um or freq>2.4653985e+14Hz
        IGM_bool = (self.wavelength<1216)
        self.f_flux[IGM_bool] = 9.1201084e-60
        self.w_flux[IGM_bool] = (9.1201082e-50 * self.cSpeed
                                /self.wavelength[IGM_bool]**2.0)

    def describeSpectrum(self):
        descriptionStr=("Fake Spectrum :" +
                        "with power law Idx = "+str(self.nuIndex)+", and "+
                        "continuum normalised to M="+
                        str(self.normalisationMag)+" at "+
                        str(self.normalisationWavelength)+" Angstroms.\n")
        emLineStr = ("Added emission lines are shown in format..\n"+
                    "'LineName':(lineCentre,intrinsicEW) : \n"+
                    str(self.emissionLineRecord)+"\n")
        detailsStr = ("Lya implemented as "+str(self.LyaType)+
                     ", and IGM absorption set to "+
                     str(self.IGM_absorptionBool)+".")
        print(descriptionStr+emLineStr+detailsStr)

    @staticmethod
    def equivWidthIntegral(lineFlux,continuumFlux,lambdaSpacing):
        """Returns EW as numerical approx. of integral.
        Given EW= integral (totalFlux-continuumFlux)/continuum 
        over wavelength, & totalFlux = lineFlux+continuumFlux.
        """
        return np.sum(lineFlux/continuumFlux) * lambdaSpacing

    @staticmethod
    def gaussianFunc(lineCenter,lineWidth,lineConst,wavelength):
        #A=norm constant, sigma&mu of gaussian evaluated @ x
        norm = lineConst / ( (2.0*np.pi)**(0.5) * lineWidth )
        exponent = -(wavelength-lineCenter)**(2.0) / (2.0 * lineWidth**2.0 )
        return norm*np.exp(exponent)
    




def main():
    """Demonstating usage of the FakeSpectrum Class"""
    nuIndex = 0.0

    baseSED = BaseSpectrum(nuIndex,100,30000,25)

    testingBool = (baseSED.wavelength>1140)&(baseSED.wavelength<1260)

    fakeSED = FakeSpectrum(nuIndex,100,30000,25,True,None)
    _ = fakeSED.showSpectrum()

    fakeSED.addEmissionLine()
    fakeSED.addEmissionLine(30, 4960.295,"OIII-I")
    fakeSED.addEmissionLine(90,5008.24,"OIII-II")   
    _ = fakeSED.showSpectrum()

    fakeSED.removeEmissionLine('lya')
    _ = fakeSED.showSpectrum()
    
    fakeSED.describeSpectrum()

if __name__ == '__main__':
    main()
