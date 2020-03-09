import numpy as np
import matplotlib.pyplot as plt

class TransmissionCurve(object):

    def __init__(self, inputFile):

        self.transmissionFile = inputFile
        self.data = self.extractContents(inputFile)

    def extractContents(self,inputFile):
        """
        Assuming ascii file format from SVO
        """
        return np.loadtxt(inputFile)


def main():
    inputFile = "TransmissionCurveFiles/Subaru_HSC.g_filter.dat"
    transCurve = TransmissionCurve(inputFile)
    print(transCurve.transmissionFile)
    print(transCurve.data)

if __name__ == '__main__':
    main()
