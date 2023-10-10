import numpy as np

from typing import Callable, Literal
from Spectrum import Spectrum
from Operator import Operator


class FunctionalOperator(Operator):
    """
    Wrapping Oprerator for any single-input function.
    """
    def __init__(self, func: Callable, label: str = None):
        self._function = func
        super().__init__(1, label)

    def __run__(self, spectra: list[Spectrum]) -> Spectrum:
        transformed = spectra[0].copy()
        transformed = self._function(spectra[0])
        return transformed


class Multiplier(Operator):

    def __init__(self):
        super().__init__(None, 'Multiplier')

    def __run__(self, spectra: list[Spectrum], *args, **kargs) -> Spectrum:
        multiplied = spectra[0].copy()
        for spectrum in spectra[1:]:
            multiplied *= spectrum
        return multiplied


class Stripper(Operator):
    """
    Strip the second spectrum from first spsectrum, and limit the minimun to be zero.
    """
    def __init__(self, min: int = 1):
        self._min = min
        super().__init__(2, 'Stripper')

    def __run__(self, spectra: list[Spectrum], *args, **kargs) -> Spectrum:
        stripped = spectra[0].copy()
        stripped = np.maximum(spectra[0] - spectra[1], self._min)
        return stripped


class Translator(Operator):
    """
    Wrapper of Numpy.roll().
    """
    def __init__(self, shift: int, reverse: bool = False):
        self._shift = shift
        self._reverse = reverse
        super().__init__(1, 'Translator')

    def __run__(self, spectra: list[Spectrum], *args, **kargs) -> Spectrum:
        translated = spectra[0].copy()
        if self._reverse:
            translated.counts = np.roll(spectra[0], -self._shift)
        else:
            translated.counts = np.roll(spectra[0], self._shift)
        return translated


class Slicer(Operator):
    """
    Wrapper of Slice operation.
    """
    def __init__(self, start, end):
        self._start = start
        self._end = end
        super().__init__(1, 'Slicer')

    def __run__(self, spectra: list[Spectrum], *args, **kargs) -> Spectrum:
        sliced = spectra[0].copy()
        sliced.counts = spectra[0][self._start: self._end]
        return sliced


class Padder(Operator):
    """
    Wrapper of Numpy.pad().
    """
    def __init__(self, length: tuple, mode: Literal['edge', 'linear_ramp', 'maximum', 'mean', 'median', 'minimum', 'reflect', 'sysmetric'] = None):
        self._length = length
        if mode is None:
            mode = 'edge'
        self._mode = mode
        super().__init__(1, 'Padder')

    def __run__(self, spectra: list[Spectrum], *args, **kargs) -> Spectrum:
        padded = spectra[0].copy()
        padded.counts = np.pad(spectra[0].counts, (self._length, self._length),
                               mode='edge')
        padded.label += self._label
        return padded
