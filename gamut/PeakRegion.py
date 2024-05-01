import numpy as np

from copy import deepcopy
from typing import Literal, Tuple
from scipy.optimize import curve_fit

from .globals import CALIBRATION_METHOD

class Calibration:

    def __init__(self, method: CALIBRATION_METHOD = 'same',
                 data: list[list[float]] | np.ndarray | None = None, params: list | None = None, label: str | None = None):
        self._method: CALIBRATION_METHOD = method

        if method == 'same':
            self._params = [0, 1]
            self._data = np.array([])
            self._by = 'params'
        elif params is not None:
            self._params = params
            self._data = np.array(data)
            self._by = 'params'
        elif isinstance(data, np.ndarray | list):
            self._data = np.array(data)
            self._params, self._fitness = self._fitparams(self._data)
            self._by = 'data'
        else:
            raise ValueError(f'Either data or params must be provided when method is {self._method}')

        if label is None:
            label = 'Calibration'
        self.label = label

    def __call__(self, inp: float | np.ndarray) -> float | np.ndarray:
        if self._method == 'same':
            return self._linear(inp, *self.params)
        elif self._method == 'linear':
            return self._linear(inp, *self.params)
        elif self._method == 'quadratic':
            return self._quadratic(inp, *self.params)
        elif self._method == 'FWHMcal':
            return self._FWHMcal(inp, *self.params)
        else:
            raise ValueError(f'Unknown curvefit method: {self._method}')

    @property
    def by(self):
        return self._by

    @property
    def method(self) -> CALIBRATION_METHOD:
        return self._method

    @property
    def params(self) -> list:
        return self._params

    @property
    def data(self) -> np.ndarray:
        return self._data

    @property
    def data_indexes(self) -> np.ndarray:
        return self._data[:, 0]

    @property
    def data_values(self) -> np.ndarray:
        return self._data[:, 1]

    @data.setter
    def data(self, data):
        self._data = data

    @property
    def label(self) -> str:
        return self._label

    @label.setter
    def label(self, label: str):
        self._label = label

    @property
    def fitness(self):
        return self._fitparams(self._data)[1]

    def width(self, inp: float | np.ndarray) -> float | np.ndarray:
        return self(inp) - self(inp-1)

    @staticmethod
    def _linear(inp: float | np.ndarray, a0: float, a1: float) -> float | np.ndarray:
        return a1 * inp + a0

    @staticmethod
    def _quadratic(inp: float | np.ndarray, a0: float, a1: float, a2: float) -> float | np.ndarray:
        return a2 * inp ** 2 + a1 * inp + a0

    @staticmethod
    def _FWHMcal(inp: float | np.ndarray, a: float, b: float, c: float) -> float | np.ndarray:
        return a + b * (inp + c * inp**2) ** 0.5

    def _fitparams(self, data: np.ndarray) -> Tuple[list[float], float]:

        if self._method == 'linear':
            func = self._linear
        elif self._method == 'quadratic':
            func = self._quadratic
        elif self._method == 'FWHMcal':
            func = self._FWHMcal
        else:
            raise ValueError(f'Unknown curvefit method: {self._method}')

        params, _ = curve_fit(func, xdata=data[:, 0], ydata=data[:, 1])
        fitted = func(self._data[:, 0], *params)
        fitness = 1 - ((fitted - self._data[:, 1])**2).sum() / ((data[:, 1] - data[:, 1].mean())**2).sum()
        return params, fitness


class Peak:
    
    def __init__(self, location: float):
        self._location = location

    @property
    def location(self):
        return self._location


class Region:

    def __init__(self, left: int, right: int, peaks: list[dict] | None = None):
        if peaks is None:
            peaks = []
        self._left = left
        self._right = right
        self.peaks = peaks
        self.slope: float = 0
        self.offset: float = 0
        self._fitted = False
        self.covariance: np.ndarray | None = None

    def __repr__(self):
        return f'Region({self.left}~{self.right})[{self.npeaks}]'        

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
    def npeaks(self):
        return len(self.peaks)

    @property
    def indexes(self):
        return np.arange(self.left, self.right+1)

    @property
    def length(self):
        return self.right - self.left + 1

    def __getitem__(self, index: int):
        return self.peaks[index]

    @staticmethod
    def _gaussian(indexes, mean, std):
        return np.exp(-(indexes - mean) ** 2 / (2 * std ** 2))

    def copy(self):
        return deepcopy(self)

    def split_peaks(self):
        splits = np.array([int((self.peaks[i+1]['location'] + self.peaks[i]['location'])/2) for i in range(len(self.peaks)-1)])
        return splits

    def fit_peak(self, peak: dict, indexes: np.ndarray | None = None) -> np.ndarray:
        if indexes is None:
            indexes = self.indexes
        fitted_peak = self._gaussian(indexes, peak['center'], peak['stderr']) * peak['height']
        return fitted_peak

    def fit_peaks(self, indexes: np.ndarray | None = None) -> np.ndarray:
        if indexes is None:
            indexes = self.indexes
        return np.array([self.fit_peak(peak, indexes) for peak in self.peaks if 'stderr' in peak.keys()]).sum(axis=0)

    def fit_baseline(self, indexes: np.ndarray | None = None) -> np.ndarray:
        if indexes is None:
            indexes = self.indexes
        fitted_baseline = indexes * self.slope + self.offset
        return fitted_baseline


# if __name__ == "__main__":
#     from Spectrum import SimulatedSpectrum, Spectrum

#     # test boundary modification
#     spec = SimulatedSpectrum()
#     peak = Region(10, 20)
#     # print(spec[peak.indexes], peak.left, peak.right)
#     # print(spec[peak.indexes], peak.left, peak.right)

#     spec.regions = []
#     spec.regions.append(peak)

#     # print([obj for obj in gc.get_referrers(peak) if isinstance(obj, Spectrum)])
#     print(peak.spectrum, peak.median('highest'), peak.counts)
