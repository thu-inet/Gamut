from Operator import Operator, DualOperator
import numpy as np


class SNIPstripper(Operator):

    def __init__(self, order):
        self._order = order
        self._label = f'-SNIPstripper[O{order}]'

    def __run__(self, spectrum):
        stripped = spectrum.copy()
        loged = np.log(np.log((spectrum.counts+1)**0.5+1)+1)
        padded = np.pad(loged, pad_width=self._order, mode='reflect')
        for p in range(self._order, 0, -1):
            compressed = np.pad(padded, pad_width=(p, 0), mode='edge')[: -p] / 2 + \
                                    np.pad(padded, pad_width=(0, p), mode='edge')[ p: ] / 2 
            padded = np.minimum(padded, compressed)
            
        recovered = (np.exp(np.exp(padded) - 1) - 1)**2
        stripped.counts -= recovered[self._order: -self._order]
        return stripped

