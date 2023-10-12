import numpy as np
from copy import deepcopy
from typing import Literal
import gc


class Calibration:

    def __init__(self, data: list, label: str = None):

        if label is None:
            label = 'Calibration'
        self.label = label

        self._data = data

    def __call__(self, inp):
        if self._data[0] is None:
            return inp
        if self._data[0] == 'linear':
            return self._linear(inp)
        elif self._data[0] == 'quadratic':
            return self._quadratic(inp)

    def _linear(self, inp):
        return self._data[1] * inp + self._data[2]

    def _quadratic(self, inp):
        return self._data[1] * inp ** 2 + self._data[2] * inp + self._data[3]

    def width(self, inp):
        if self._data[0] is None:
            return 1
        if self._data[0] == 'linear':
            return self._data[1]


class Region:

    def __init__(self, left: int, right: int):

        self._left = left
        self._right = right

    @property
    def left(self):
        return self._left

    @left.setter
    def left(self, left: int):
        self._left = left

    @property
    def right(self):
        return self._right

    @right.setter
    def right(self, right: int):
        self._right = right

    @property
    def indexes(self):
        return np.arange(self.left, self.right+1)

    @property
    def length(self):
        return self.right - self.left + 1

    @property
    def spectrum(self):
        # if not hasattr(self, '_spectrum'):
        from Spectrum import Spectrum
        _spectrum = next((obj for obj in gc.get_objects() if isinstance(obj, Spectrum) and hasattr(obj, 'peaks') and (self in obj.peaks)), None)
        return _spectrum

    @property
    def counts(self):
        return self.spectrum.counts[self.indexes]

    def median(self, mode: Literal["center", "highest"] = 'center'):
        if mode == 'center':
            median = round((self.left + self.right) / 2)
        elif mode == 'highest':
            median = np.argmax(self.counts) + self.left
        return median

    def copy(self):
        return deepcopy(self)


if __name__ == "__main__":
    from Spectrum import SimulatedSpectrum, Spectrum

    # test boundary modification
    spec = SimulatedSpectrum()
    peak = Region(10, 20)
    # print(spec[peak.indexes], peak.left, peak.right)
    # print(spec[peak.indexes], peak.left, peak.right)

    spec.peaks = []
    spec.peaks.append(peak)

    # print([obj for obj in gc.get_referrers(peak) if isinstance(obj, Spectrum)])
    print(peak.spectrum, peak.median('highest'), peak.counts)