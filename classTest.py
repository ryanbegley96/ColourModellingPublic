class BaseClass(object):
    
    def __init__(self,a,b):
        self.attribute1 = a
        self.attribute2 = b

    def showBaseAttributes(self):
        print("Attribute1 = {}, Attribute2 = {}".format(self.attribute1, 
                                                        self.attribute2))

class SubClass(BaseClass):

    def __init__(self,a=1,b=10,c=100):
        self.attribute3 = (c,"helloworld!")
        self.attribute4 = self.attribute3**2.0
        super().__init__(a,b)

    @property
    def attribute3(self):
        return self._attribute3 
    
    @attribute3.setter
    def attribute3(self,valTuple):
        """
        note - arguement in setter takes 1 value only so pass a tuple, but
        want to test way of setting default to none for case of flux changes as 
        need to call fluxdensityUnitSwap with a wavelength but give option of wl.
        Then would assign using a tuple.
        """

        val,test = valTuple
        self._attribute3 = val
        self._attribute4 = val**2.0
        print("Testing use of test tuple in attribute3.setter : "+test)

    @property
    def attribute4(self):
        return self._attribute4

    @attribute4.setter
    def attribute4(self,val):
        self._attribute3 = val**0.5    
        self._attribute4 = val

    def showSubAttributes(self):
        print("Attribute3 = {}, Attribute4 = {}".format(self.attribute3, 
                                                        self.attribute4))

    def __iadd__(self,val):
        self.attribute3+= val[0]                                                
        return self

class SomeClass(object):
    def __init__(self, n):
        self.list = range(0, n)

    @property
    def list(self):
        return self._list
    @list.setter
    def list(self, val):
        self._list = val
        self._listsquare = [x**2 for x in self._list ]

    @property
    def listsquare(self):
        return self._listsquare
    @listsquare.setter
    def listsquare(self, val):
        self.list = [int(pow(x, 0.5)) for x in val]


def main():
    a,b,c = 0.0,1.0,2.0
    baseClass = BaseClass(a,b)
    baseClass.showBaseAttributes()

    subClass = SubClass(a,b,c)
    subClass.showSubAttributes()
    print("\n")

    testStr = "HelloWorld!"
    print("New test string to be :"+testStr)
    subClass.attribute3 = (10,testStr)
    print("\nNow showing the new attribute 3:")
    subClass.showSubAttributes()
    
    subClass.attribute4 = 1000
    subClass.showSubAttributes()

if __name__ == '__main__':
    main()

# import copy 
# import numpy as np
# import matplotlib.pyplot as plt

# class BaseSpectrum(object):

#     def __init__(self,nuIndex,lowerLambda,upperLambda,normalisationMag):
        
#         self.cSpeed = 3.0E+8
#         self.normalisationWavelength = 25000.0 #2.5microns
#         self.lambdaSpacing = 1

#         self.nuIndex = nuIndex
#         self.nuNorm = self.setNuNorm(normalisationMag)

#         self.lambdaIndex = nuIndex-2.0
#         self.lambdaNorm = self.setLambdaNorm(normalisationMag)

#         self.wavelength = np.arange(lowerLambda,upperLambda,self.lambdaSpacing)
#         self.frequency = 1.0E-10 * self.cSpeed /self.wavelength

#         self.f_flux = self.spectrumNuModel(self.frequency,self.nuIndex)
#         self.w_flux = self.spectrumLambdaModel(self.wavelength,
#                                                self.lambdaIndex)

#     # @property
#     # def f_flux(self):
#     #     return self._f_flux
#     # @f_flux.setter
#     # def f_flux(self,newFlux):
#     #     self._f_flux = newFlux
#     #     self._w_flux = 


#     def setLambdaNorm(self,normalisationMag):
#         w_flux = self.fluxdensityUnitSwap( self.MagToFlux(normalisationMag),
#                                     'frequency', self.normalisationWavelength )
#         return w_flux / self.normalisationWavelength**self.lambdaIndex

#     def setNuNorm(self,normalisationMag): #3.6307805E-30
#         self.normalisationNu = 1E10*self.cSpeed/self.normalisationWavelength
#         return (self.MagToFlux(normalisationMag)
#                 / self.normalisationNu**self.nuIndex )

#     def fluxdensityUnitSwap(self,fluxdensity,currentUnit='wavelength',
#                         wavelength=None):
#         """
#         Converts from F_lam [erg/s/cm2/A] <-> F_nu [erg/s/cm2/Hz]
#         """
#         if wavelength is None:
#             wavelength = self.wavelength

#         if currentUnit == 'wavelength':  #fv=3.34x10^-19 * l**2 * fl
#             return fluxdensity * wavelength**2.0 * 1.0E-10 / self.cSpeed
#         elif currentUnit == 'frequency': #fl=3x10^18 * fv / l**2
#             return fluxdensity * 1.0E+10 * self.cSpeed / wavelength**2.0
#         else:
#             raise ValueError( ("Bad Argument in fluxdensityUnitSwap, currentUnit must"
#                             "be 'wavelength' or 'frequency'.")
#                         )

#     def spectrumNuModel(self,frequency,nuIndex):
#         """
#         Spectrum modelled as fv=a*frequency**nuIndex
#         """
#         return self.nuNorm*frequency**nuIndex

#     def spectrumLambdaModel(self,wavelength,lambdaIndex):
#         """
#         Spectrum modelled as fl=b*wavelength**alpha, where alpha=nuIndex-2
#         """
#         return self.lambdaNorm*wavelength**lambdaIndex

#     def f_fluxToAB(self,f_flux):
#         return -2.5*np.log10(f_flux)-48.6
    
#     def w_fluxToAB(self,w_flux,wavelength=None):
#         if wavelength is None:
#             return -2.5*np.log10(self.fluxdensityUnitSwap(w_flux))-48.6
#         else:
#             return -2.5*np.log10(self.fluxdensityUnitSwap(w_flux,'wavelength',
#                                 wavelength) )-48.6

#     def MagToFlux(self,mag):
#         return 10.0**(-(mag+48.6)/2.5)

#     def showSpectrum(self):
#         _,axs = plt.subplots()
#         axs.plot(self.wavelength/1E4,self.f_fluxToAB(self.f_flux))
#         axs.set_xlabel(r'$\lambda/\mu$m')
#         axs.set_ylabel(r'$m_{AB}$')
#         axs.set_xlim(xmin=0.,xmax=3.0)
#         axs.set_ylim(ymax=30,ymin=22)
#         axs.invert_yaxis()
#         # plt.savefig("fakeSpectrum.png")
#         plt.show()

#     def testMessage(self):
#         print("Done BaseClass.")   

# class FakeSpectrum(BaseSpectrum):

#     def __init__(self,nuIndex=0, lowerLambda=100,upperLambda=30000,
#                  normalisationMag=25,activateIGM=True,activateLya=None):
#         super().__init__(nuIndex,lowerLambda,upperLambda,normalisationMag)
#         if activateLya is "Dirac":
#             self.implementDiracLya()    
#         if activateLya is "Gaussian":
#             self.implementGaussianLya()

#         if activateIGM:
#             self.IGM_absorption()

#      ### Testing Region  ->  Implementing fake emission line
#     ########################################################################
#     @staticmethod
#     def gaussianFunc(lineCenter,lineWidth,lineConst,wavelength):
#         #A=norm constant, sigma&mu of gaussian evaluated @ x
#         norm = lineConst / ( (2.0*np.pi)**(0.5) * lineWidth )
#         exponent = -(wavelength-lineCenter)**(2.0) / (2.0 * lineWidth**2.0 )
#         return norm*np.exp(exponent)
    
#     def implementGaussianLya(self):
#         equivalentWidth = 400.0
#         lineCenter = 1215.7 #wavelength of Lya emission line
#         lineWidth = 2.0 #sigma of Gaussian in Angstroms
#         integLimits = (lineCenter-5.*lineWidth,lineCenter+5.*lineWidth)
        
#         idx = np.argmin(np.abs(self.wavelength-1215.7)) #fc val @ 1216A
#         lineConst = equivalentWidth * self.w_flux[idx] #7.366E-18#erg/s/cm2/A

#         wavelengthBool = ( (self.wavelength > integLimits[0]) & 
#                             (self.wavelength < integLimits[1]) )
#         lineWavelengths = self.wavelength[wavelengthBool]
#         old_w_flux = copy.deepcopy(self.w_flux)[wavelengthBool]

#         lineFlux = self.gaussianFunc(lineCenter,lineWidth,lineConst,
#                                      lineWavelengths)
#         self.w_flux[wavelengthBool] += lineFlux
#         self.f_flux[wavelengthBool] += self.fluxdensityUnitSwap(lineFlux,
#                                             'wavelength',lineWavelengths)
#         equivalentWidthCalc = np.sum( lineFlux/old_w_flux ) * self.lambdaSpacing
#         print("Equivalent Width :",equivalentWidthCalc,"/Angstroms")
        

#     def implementDiracLya(self):
#         idx = np.argmin(np.abs(self.wavelength-1215.7))
#         equivalentWidth = 200.0
#         A = equivalentWidth * self.w_flux[idx] #7.366E-18
#         print("Old w flux :",self.w_flux[idx])
#         print("Old f flux :",self.f_flux[idx])
#         self.w_flux[idx] += A
#         print("New w flux :",self.w_flux[idx])
#         print("New f flux (after w flux update):",self.f_flux[idx])
#         self.f_flux[idx] += self.fluxdensityUnitSwap(A,'wavelength',1216.0)
#         print("New f flux (after f flux update):",self.f_flux[idx])  
#     ########################################################################


#     def IGM_absorption(self):
#         #Flux should be set to 0 at lambda<1.216um or freq>2.4653985e+14Hz
#         IGM_bool = (self.wavelength<1216)
#         self.f_flux[IGM_bool] = 9.1201084e-60
#         self.w_flux[IGM_bool] = (9.1201082e-50 * self.cSpeed
#                                 /self.wavelength[IGM_bool]**2.0)
#     def testMessage(self):
#         print("Done SubClass.")   


# def main():
#     nuIndex = 0.0

#     baseSED = BaseSpectrum(nuIndex,100,30000,25)
#     baseSED.testMessage()

#     fakeSED = FakeSpectrum(nuIndex,100,30000,25,True,'Gaussian')
#     # fakeSED.showSpectrum()

# if __name__ == '__main__':
#     main()
