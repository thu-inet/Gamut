import numpy as np
import matplotlib.pyplot as plt

from pathlib import Path
from copy import deepcopy
from typing import Callable, Literal

import plot_settings
from PeakRegion import PeakRegion


class Spectrum(np.ndarray):
    """
    Spectrum data.
    """
    def __new__(cls, counts: np.ndarray, *args, **kargs):
        self = np.asarray(counts).view(cls)
        return self

    def __array_finalize__(self, obj):
        if obj is None:
            return
        self._label = getattr(obj, 'label', None)
        self._peaks = getattr(obj, 'peaks', None)
        if (attrs := getattr(obj, 'attrs', None)) is not None:
            self._attrs = deepcopy(attrs)
        else:
            self._attrs = None

    def __init__(self, counts: list | np.ndarray, label: str = None, attrs: dict = None):
        if label is None:
            label = 'Spectrum'
        self._label = label

        if attrs is None:
            attrs = {}
        self._attrs = attrs

    def __array_wrap__(self, out_arr, context=None):
        if out_arr.ndim == 0:
            return out_arr.item()
        else:
            return np.ndarray.__array_wrap__(self, out_arr, context)

    def __repr__(self):
        return f"Spectrum[{self.label}]\n \
|Shape: {self.shape}\n \
|Head:  {list(self[:min(5, len(self))])}\n \
|Attrs: {self.attrs.keys()}\n"

    def __str__(self):
        return f"Spectrum[{self.label}]"

    @property
    def label(self):
        return self._label

    @label.setter
    def label(self, value):
        self._label = value

    @property
    def indexes(self):
        return np.arange(len(self))

    @property
    def attrs(self):
        return self._attrs

    @attrs.setter
    def attrs(self, value):
        self._attrs = value

    @property
    def counts(self):
        return np.asarray(self)

    @property
    def length(self):
        return self.shape[0]

    @property
    def peaks(self):
        return self._peaks

    @peaks.setter
    def peaks(self, peaks: list[PeakRegion]):
        self._peaks = peaks

    def copy(self):
        return deepcopy(self)

    def area_estimate(self, peak: PeakRegion) -> tuple:
        """
        Estimate peak area of spectrum
        """
        total = sum(self[peak.indexes])
        baseline = (self[peak.left] + self[peak.right]) * peak.length / 2
        return total, baseline

    def plot(self, *args, axes: plt.Axes = None, **kargs) -> plt.Axes:
        """
        Plot the spectrum.
        """
        if axes is None:
            axes = plt.gca()
        if 'label' not in kargs:
            kargs['label'] = self.label
        axes.plot(self, *args, **kargs)
        axes.set_ylim(0, )
        axes.set_xlim(0, )
        return axes

    def plot_peaks(self, *args, axes: plt.Axes = None, **kargs) -> plt.Axes:
        if not hasattr(self, 'peaks'):
            raise ValueError(f"{self} does not have peaks.")
        if axes is None:
            axes = plt.gca()
        for i, peak in enumerate(self.peaks):
            axes.plot(peak.indexes, self[peak.indexes], *args, **kargs, marker=plot_settings.markers(i))
    
    # def plot_fitted_peaks(self, *args, axes: plt.Axes = None, **kargs) -> plt.Axes:
    #     if not hasattr(self, 'peaks'):
    #         raise ValueError(f"{self} does not have peaks.")
    #     if axes is None:
    #         axes = plt.gca()
    #     for i, peak in enumerate(self.peaks):
    #         if np.any([not hasattr(peak, attr) for attr in ['params', 'shapes', 'heights', 'fcounts']]):
    #             continue
    #         axes.plot(peak.indexes, peak.fcounts, *args, **kargs, marker=plot_settings.markers)
    #         axes.fill_between()

    @classmethod
    def from_GammaVision(cls, filename: str):
        if Path(filename).suffix.lower() not in ['.spe', '.chn']:
            raise ValueError(f"{filename} is not a valid GammaVision file.")
        with open(filename, 'r') as fopen:
            filelines = fopen.readlines()
        data_index = next((i for i, line in enumerate(filelines) if line.startswith('$DATA')), None)
        index = data_index + 2
        counts = []
        while filelines[index].strip().isdigit():
            counts.append(int(filelines[index].strip()))
            index += 1
        self = cls(counts)
        return self
        

class SimulatedSpectrum(Spectrum):
    """
    Simulated Spectrum Generator for algorithm test.
    """
    def __new__(cls, length: int = 200, *args, **kargs):
        counts = np.zeros(length)
        self = super().__new__(cls, counts)
        return self

    def __init__(self,
                 length: int = 200,
                 peaks_info: list = [[30, 2.5, 4000], [70, 3, 2000], [110, 4, 1000],  [150, 4.5, 2400], [165, 4.3, 1600]],
                 base_intensity: int = 100,
                 base_amplitude: int = 100,
                 base_function: Callable = lambda x: x,
                 label: str = None):

        self.peaks_info = peaks_info
        self.base_intensity = base_intensity
        self.base_amplitude = base_amplitude
        self.base_function = base_function

        peaks = []
        for peak_info in peaks_info:
            peaks.append(PeakRegion(int(peak_info[0]-peak_info[1]*3), int(peak_info[0]+peak_info[1]*3)))

        self.baseline = self._gen_base()
        self.peakform = self._gen_peaks()
        self[:] = self.baseline[:] + self.peakform[:]
        super().__init__(self.baseline + self.peakform, label)

    def _gen_base(self):
        base = self.base_function(np.linspace(0, 1, self.length))
        base = (base - base.min()) / (base.ptp() + 1E-3)
        base = base * self.base_amplitude + self.base_intensity
        base = np.random.normal(loc=base, scale=base**0.5, size=(self.length,))
        return base

    def _gen_peaks(self):
        peak = np.zeros((self.length,))
        for peak_info in self.peaks_info:
            peak += self._gen_peak(*peak_info)
        return peak

    def _gen_peak(self, centroid, stderror, area):
        channels = np.arange(0, self.length)
        amplitude = area / stderror / (2*np.pi)**0.5
        return amplitude * np.exp(- (channels-centroid)**2 / stderror**2 / 2)


simuspecs = {
    'simple': SimulatedSpectrum(base_intensity=100, base_amplitude=100, base_function=lambda x:  x**0 + 1,
                                    peaks_info=[[95, 6, 4000]]),
    'doublepeak_slight': SimulatedSpectrum(base_intensity=100, base_amplitude=100, base_function=lambda x:  x,
                                    peaks_info=[[70, 6, 4000], [105, 7, 3000]]),
    'doublepeak_normal_narrow': SimulatedSpectrum(base_intensity=100, base_amplitude=1000, base_function=lambda x:  1 / (x + 0.1),
                                    peaks_info=[[90, 3, 4000], [100, 2.5, 3000]]),
    'doublepeak_normal': SimulatedSpectrum(base_intensity=100, base_amplitude=100, base_function=lambda x:  1 / (x + 1),
                                    peaks_info=[[70, 6, 4000], [100, 10, 3000]]),
    'doublepeak_severe': SimulatedSpectrum(base_intensity=200, base_amplitude=100, base_function=lambda x:  np.abs((x-0.4)),
                                    peaks_info=[[70, 7, 4000], [90, 10, 3000]]),
    'synthesized': SimulatedSpectrum(base_intensity=100, base_amplitude=100, base_function=lambda x:  x,
                                    peaks_info=[[30, 2.5, 4000], [70, 3, 2000], [110, 4, 1000],  [150, 4.5, 2400], [165, 4.3, 1600]])
}


if __name__ == "__main__":

    # Initilization
    spec = Spectrum(np.arange(1, 8), label="test", attrs={"a": 1, "b": 2})

    # Check behavior
    print(f"1. Print length = {len(spec)}")
    print(f"2. Slice like a array: spec[1:3] = {spec[1:3]}, spec[1:] = {spec[1:]}, spec[1:3:2] = {spec[1:3:2]}")
    print(f"3. Calculate like a array: spec+1 = {(spec+1)}, spec+spec = {(spec+spec)}, sum(spec) = {sum(spec)}")

    # test spectrum modification
    spec = Spectrum(np.zeros(10))
    print(id(spec), spec.sum(), spec.counts)
    spec[:5] = 1
    print(id(spec), spec.sum(), spec.counts)
    spec[5:] = 2
    print(id(spec), spec.sum(), spec.counts)
    spec[:] = np.arange(10)
    print(id(spec), spec.sum(), spec.counts)

    # test spectrum
    simu = SimulatedSpectrum()
    simu.plot()
    plt.savefig("fig.jpg")
