from Operator import Operator, DualOperator
import numpy as np


class Translator(Operator):

    def __init__(self, shift, reverse=False):
        self._shift = shift
        self._reverse = reverse
        self._label = f'-shift[{self._shift}]'

    def __run__(self, spectrum):
        trslted = spectrum.copy()
        if self._reverse:
            trslted.counts = np.roll(spectrum, -self._shift)
        else:
            trslted.counts = np.roll(spectrum, self._shift)
        trslted.label += self._label
        return trslted


class Slicer(Operator):

    def __init__(self, start, end):
        self._start = start
        self._end = end
        self._label = f'-slice[{self._start}:{self._end}]'

    def __run__(self, spectrum):
        slced = spectrum.copy()
        slced.counts = spectrum[self._start: self._end]
        slced.label += self._label
        return slced


class Extender(Operator):

    def __init__(self, length):
        self._length = length
        self._label = f'-extend[{self._length}]'

    def __run__(self, spectrum):
        extded = spectrum.copy()
        extded.counts = np.pad(spectrum.counts, (self._length, self._length),
                               mode='edge')
        extded.label += self._label
        return extded

