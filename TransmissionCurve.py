import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

class TransmissionCurve(object):

    def __init__(self, inputFile):

        self.transmissionFile = inputFile
        self.wavelength,self.transmission = self.extractContents(inputFile)

    def extractContents(self,inputFile):
        """
        Assuming ascii file format from SVO
        """
        data = np.loadtxt(inputFile)
        return data[:,0],data[:,1]

    def returnInterpolationFunc(self):
        """
        Creates interpolation function which can be used to sample.
        """
        self.interpFunc = interp1d(self.wavelength,self.transmission,
                                    kind='cubic')
        return self.interpFunc

    def showTransmissionCurve(self):
        _,axs = plt.subplots()
        axs.plot(self.wavelength/1E4,self.transmission)
        axs.set_xlabel(r"$\lambda/\mu$m")
        axs.set_ylabel("Transmission/%")
        axs.set_ylim(ymin=0)
        return axs

def main():
    inputFile = "TransmissionCurveFiles/Subaru_HSC.Y_filter.dat"
    transCurve = TransmissionCurve(inputFile)

    interpFunc = transCurve.returnInterpolationFunc()
    xData = np.linspace(np.min(transCurve.wavelength),
                        np.max(transCurve.wavelength),1000,endpoint=True)
    yData = interpFunc(xData)
    axs = transCurve.showTransmissionCurve()
    axs.plot(xData/1E4,yData)
    plt.show()

if __name__ == '__main__':
    main()





"""
Experimental Code:
>Needs high order > 100 to fit decently but this is time intensive.

from symfit import parameters,variables,sin,cos,Fit

def fourier_series(x, f, n=0):

    #Returns a symbolic fourier series of order `n`.

    #:param n: Order of the fourier series.
    #:param x: Independent variable
    #:param f: Frequency of the fourier series

    # Make the parameter objects for all the terms
    a0, *cos_a = parameters(','.join(['a{}'.format(i) for i in range(0, n + 1)]))
    sin_b = parameters(','.join(['b{}'.format(i) for i in range(1, n + 1)]))
    # Construct the series
    series = a0 + sum(ai * cos(i * f * x) + bi * sin(i * f * x)
    for i, (ai, bi) in enumerate(zip(cos_a, sin_b), start=1))
    return series

x, y = variables('x, y')
w, = parameters('w')
model_dict = {y: fourier_series(x, f=w, n=125)}
print(model_dict)
fit = Fit(model_dict, x=transCurve.wavelength/1E4, y=transCurve.transmission)
fit_result = fit.execute()
print(fit_result)

_,axs = plt.subplots()
axs.plot(transCurve.wavelength/1E4,transCurve.transmission)
axs.plot(transCurve.wavelength/1E4, fit.model(x=transCurve.wavelength/1E4, **fit_result.params).y,ls=':')
axs.set_xlabel(r"$\lambda/\mu$m")
axs.set_ylabel("Transmission/%")
axs.set_ylim(ymin=0)
plt.show()

"""
