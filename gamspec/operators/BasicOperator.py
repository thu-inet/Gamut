'''
Author: albertzhang albert.zhangweij@outlook.com
Date: 2023-12-25 10:20:01
Description: 

Copyright (c) 2023 by THU-RSAG, All Rights Reserved. 
'''
import numpy as np

from typing import Callable, Literal
from ..Spectrum import Spectrum
from ..Operator import Operator


class Passer(Operator):

    def __init__(self):
        super().__init__(1, 'Passer')

    def __run__(self, spectra: list[Spectrum], *args, **kargs) -> Spectrum:
        return spectra[0].copy()


class Stripper(Operator):
    """
    Strip the second spectrum from first spsectrum, and limit the minimun to be zero.
    """
    def __init__(self, min: int = 0):
        self._min = min
        super().__init__(2, 'Stripper')

    def __run__(self, spectra: list[Spectrum], *args, **kargs) -> Spectrum:
        stripped = spectra[0].copy()
        stripped[:] = np.maximum(spectra[0] - spectra[1], self._min)
        return stripped


class Functionor(Operator):
    """
    Wrapping Oprerator for any single-input function.
    """
    def __init__(self, func: Callable, label: str = 'Functionor'):
        self._function = func
        super().__init__(1, label)

    def __run__(self, spectra: list[Spectrum]) -> Spectrum:
        transformed = spectra[0].copy()
        transformed = self._function(spectra[0])
        return transformed


class Combinor(Operator):

    def __init__(self, inp_num: int, func: Callable, label: str = 'Combinor'):
        self.func = func
        super().__init__(inp_num, label)

    def __run__(self, spectra: list[Spectrum], *args, **kargs) -> Spectrum:
        combined = self.func(*spectra)
        return combined


class Convolutionor(Operator):
    
    def __init__(self, kernel: list[float], mode: Literal['same', 'full', 'valid'] = 'same', label: str = 'Convolutionor'):
        self._kernel = np.array(kernel)
        self._mode = mode
        super().__init__(1, label)
    
    def __run__(self, spectra: list[Spectrum], *args, **kargs) -> Spectrum:
        convolved = spectra[0].copy()
        convolved[:] = np.convolve(spectra[0], self._kernel, mode=self._mode)
        return convolved


class Multiplier(Operator):

    def __init__(self):
        super().__init__(None, 'Multiplier')

    def __run__(self, spectra: list[Spectrum], *args, **kargs) -> Spectrum:
        multiplied = spectra[0].copy()
        for spectrum in spectra[1:]:
            multiplied[:] *= spectrum
        return multiplied


class Translator(Operator):

    def __init__(self, shift: int, reverse: bool = False):
        self._shift = shift
        self._reverse = reverse
        super().__init__(1, 'Translator')

    def __run__(self, spectra: list[Spectrum], *args, **kargs) -> Spectrum:
        translated = spectra[0].copy()
        if self._reverse:
            translated[:] = np.roll(spectra[0], -self._shift)
        else:
            translated[:] = np.roll(spectra[0], self._shift)
        return translated






