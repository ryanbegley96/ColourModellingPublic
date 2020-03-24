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
        axs.plot(self.wavelength/1E4,self.transmission,label="Real")
        axs.set_xlabel(r"$\lambda/\mu$m")
        axs.set_ylabel("Transmission/%")
        axs.set_ylim(ymin=0)
        return axs

def main():
    inputFile = "TransmissionCurveFiles/Subaru_HSC.Y.dat"
    # inputFile = "TransmissionCurveFiles/Paranal_VISTA.Y.dat"
    transCurve = TransmissionCurve(inputFile)

    interpFunc = transCurve.returnInterpolationFunc()
    xData = np.linspace(np.min(transCurve.wavelength),
                        np.max(transCurve.wavelength),1000,endpoint=True)
    yData = interpFunc(xData)
    axs = transCurve.showTransmissionCurve()
    axs.plot(xData/1E4,yData,label="Interpolated")
    axs.legend(fontsize=10)
    plt.show()

if __name__ == '__main__':
    main()
