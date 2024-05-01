import numpy as np

from copy import deepcopy
from typing import Literal
from scipy.optimize import curve_fit


class Calibration:

    def __init__(self, method: Literal['same', 'linear', 'quadratic', 'FWHMcal'] = 'same',
                 data: list = None, params: list = None, label: str = None):
        self._method = method
        if method == 'same':
            self._params = [0, 1, 0]
        else:
            if data is None and params is None:
                raise ValueError(f'Either data or params must be provided when method is {self._method}')
            elif params is not None:
                self._params = params
            else:
                if len(data[0]) <= 1:
                    raise ValueError(f'Data must have at least two points when method is {self._method}')
                self._data = np.array(data)
                self._params, self._fitness = self._curvefit()

        if label is None:
            label = 'Calibration'
        self.label = label

    def __call__(self, inp: float | np.ndarray) -> float | np.ndarray:
        if self._method == 'same':
            return self._linear(inp, *self._params)
        elif self._method == 'linear':
            return self._linear(inp, *self._params)
        elif self._method == 'quadratic':
            return self._quadratic(inp, *self._params)
        elif self._method == 'FWHMcal':
            return self._FWHMcal(inp, *self._params)

    @property
    def method(self) -> Literal['same', 'linear', 'quadratic', 'FWHMcal']:
        return self._method

    @property
    def params(self) -> list:
        return self._params

    @property
    def data(self) -> list:
        return self._data

    @property
    def label(self) -> str:
        return self._label

    @label.setter
    def label(self, label: str):
        self._label = label

    @staticmethod
    def _linear(inp: float | np.ndarray, a0: float, a1: float, a2: float) -> float | np.ndarray:
        return a1 * inp + a0

    @staticmethod
    def _quadratic(inp: float | np.ndarray, a0: float, a1: float, a2: float) -> float | np.ndarray:
        return a2 * inp ** 2 + a1 * inp + a0

    @staticmethod
    def _FWHMcal(inp: float | np.ndarray, a: float, b: float, c: float) -> float | np.ndarray:
        return a + b * (inp + c * inp**2) ** 0.5

    def _curvefit(self) -> list[float]:
        if self._method == 'linear':
            func = self._linear
        elif self._method == 'quadratic':
            func = self._quadratic
        elif self._method == 'FWHMcal':
            func = self._FWHMcal
        else:
            raise ValueError(f'Unknown curvefit method: {self._method}')
        params, _ = curve_fit(func, xdata=self._data[0, :], ydata=self._data[1, :])
        fitted = func(self._data[0, :], *params)
        fitness = 1 - ((fitted - self._data[1, :])**2).sum() / ((self._data[1, :] - self.data[1, :].mean())**2).sum()
        return params, fitness

    @property
    def width(self, inp):
        return self(inp) - self(inp-1)

    @property
    def data(self):
        return self._data

    @data.setter
    def data(self, data):
        self._data = data


class Region:

    def __init__(self, left: int, right: int, peaks: list[dict] = None):
        self._left = left
        self._right = right
        if peaks is None:
            peaks = []
        self._peaks = peaks

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
    def peaks(self):
        return self._peaks

    @peaks.setter
    def peaks(self, peaks: list[dict]):
        self._peaks = peaks

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
        return self._peaks[index]

    @staticmethod
    def _gaussian(indexes, mean, std):
        return np.exp(-(indexes - mean) ** 2 / (2 * std ** 2))

    def copy(self):
        return deepcopy(self)

    def split_peaks(self):
        splits = np.array([int((self.peaks[i+1]['location'] + self.peaks[i]['location'])/2) for i in range(len(self.peaks)-1)])
        return splits

    def fit_peak(self, peak: dict, indexes: np.ndarray = None) -> np.ndarray:
        if indexes is None:
            indexes = self.indexes
        fitted_peak = self._gaussian(indexes, peak['center'], peak['stderr']) * peak['height']
        return fitted_peak

    def fit_peaks(self, indexes: np.ndarray = None) -> np.ndarray:
        if indexes is None:
            indexes = self.indexes
        return np.array([self.fit_peak(peak, indexes) for peak in self.peaks]).sum(axis=0)

    def fit_baseline(self, indexes: np.ndarray = None) -> np.ndarray:
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
