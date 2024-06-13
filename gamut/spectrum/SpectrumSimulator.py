import numpy as np
import matplotlib.pyplot as plt
import xml.etree.ElementTree as ET
import xml.dom.minidom as minidom
import pickle
import pandas as pd

from collections import namedtuple
from time import strftime, localtime
from pathlib import Path
from copy import deepcopy
from typing import Callable, Any
from numpy.typing import NDArray
from xml.dom.minidom import parse as parse_xml
from re import match


from ..classes.Calibration import Calibration
from ..globals import StructuredFunction
from .Spectrum import Spectrum


class SpectrumSimulator:
    
    def __init__(self, baseline_function: StructuredFunction,
                 peaks_function: StructuredFunction,
                 noise_function: StructuredFunction):
        self._baseline_function = baseline_function
        self._peaks_function = peaks_function
        self._noise_function = noise_function

    def generate_baseline(self, length: int, params: list[float], *args, **kargs) -> Spectrum:
        indexes = np.arange(length)
        base = self._baseline_function.function(indexes, params)
        return Spectrum(base, *args, **kargs)
    
    def generate_peaks(self, length: int, peaks_info: list[list[float]], *args, **kargs) -> Spectrum:
        indexes = np.arange(length)
        peaks = np.zeros(length)
        for peak_info in peaks_info:
            peaks += self._peaks_function.function(indexes, peak_info)
        return Spectrum(peaks, *args, **kargs)

    def add_noise(self, spectrum: Spectrum, params: list[float], rand_seed: int = 100) -> Spectrum:
        noised_spectrum = self._noise_function.function(spectrum, params)
        return noised_spectrum

    def generate(self,
                 length: int,
                 baseline_params: list[float],
                 peaks_info: list[list[float]],
                 noise_params: list[float],
                 *args, **kargs) -> Spectrum:
        baseline = self.generate_baseline(length, baseline_params, *args, **kargs)
        peaks = self.generate_peaks(length, peaks_info, *args, **kargs)
        spectrum = baseline + peaks
        noised_spectrum = self.add_noise(spectrum, noise_params)
        return noised_spectrum

spectrum_parameter = namedtuple('spectrum_parameters', ['label', 'length', 'base_params', 'peaks_params'])
spectrum_parameters = {
    'simple': spectrum_parameter('simple', 100, [100, 100], [[95, 6, 4000]]),
}

simuspecs = {
    'simple': {'base_intensity': 100, 'base_amplitude': 100, 'label': 'simple', 'length': 200,   # FWHM~2.355*5~14.13, params[linear]=[0.14873, 0]
               'base_function': lambda x:  x**0 + 1, 'peaks_info': [[95, 6, 4000]], 'FWHMcal': Calibration('linear', params=[0, 0.14873, 0])},
    'double_slight': {'base_intensity': 100, 'base_amplitude': 100, 'label': 'double_slight', 'length': 200,  
                      'base_function': lambda x:  x, 'peaks_info': [[70, 6, 4000], [105, 7, 3000]], 'FWHMcal': Calibration('linear', params=[9.4225, 0.067285, 0])},
    'double_normal_narrow': {'base_intensity': 100, 'base_amplitude': 1000, 'label': 'double_normal_narrow', 'length': 200,  
                             'base_function': lambda x:  1 / (x + 0.1), 'peaks_info': [[90, 3, 4000], [100, 2.5, 3000]]},
    'double_normal': {'base_intensity': 100, 'base_amplitude': 100, 'label': 'double_normal', 'length': 200,
                      'base_function': lambda x:  1 / (x + 1), 'peaks_info': [[70, 6, 4000], [100, 10, 3000]]},
    'double_severe': {'base_intensity': 200, 'base_amplitude': 100, 'label': 'double_severe', 'length': 200,  
                      'base_function': lambda x:  np.abs((x-0.4)), 'peaks_info': [[70, 7, 4000], [90, 10, 3000]]},
    'synthesized': {'base_intensity': 100, 'base_amplitude': 100, 'label': 'synthesized_multiplets', 'length': 200,  
                    'base_function': lambda x:  x, 'peaks_info': [[30, 2.5, 4000], [60, 3, 2000], [100, 4, 1000],  [130, 4.5, 2400], [145, 4.3, 1600], [175, 4.7, 2600]]},
    'synthesized_2': {'base_intensity': 100, 'base_amplitude': 100, 'label': 'synthesized_multiplets', 'length': 500,  # with FWHM cal: 2+0.03*(E+0.05*E**2)**0.5
                      'base_function': lambda x:  x, 'peaks_info': [[30, 2.33, 2000], [45, 2.42, 2000], [70, 2.57, 4000], [75, 2.61, 1500],
                                                                    [110, 2.83, 1500], [117, 2.87, 500], [126, 2.90, 2000],
                                                                    [150, 3.07, 2400], [165, 3.16, 1600],
                                                                    [230, 3.56, 4000], [275, 3.83, 2000], [290, 3.93, 2500], [305, 4.02, 1500],
                                                                    [310, 4.04, 1500], [317, 4.1, 500], [344, 4.22, 2000], [364, 4.40, 2000], [386, 4.50, 2000],
                                                                    [450, 4.91, 2400], [470, 5.01, 1600]]},
    'single_peaks': {'base_intensity': 100, 'base_amplitude': 100, 'label': 'synthesized_singlets', 'length': 500,  # with FWHM cal: 2+0.03*(E+0.05*E**2)**0.5
                     'base_function': lambda x:  x, 'peaks_info': [[110, 2.95, 2000], [130, 3, 2000], [180, 3.41, 4000], [230, 3.75, 1500],
                                                                   [285, 4.14, 2500], [335, 4.4, 2500], [375, 4.7, 1000],
                                                                   [410, 4.95, 2400], [450, 5.2, 1600]]},
    'single_peaks_bigger': {'base_intensity': 100, 'base_amplitude': 100, 'label': 'synthesized_singlets', 'length': 500,  # with FWHM cal: 2+0.03*(E+0.05*E**2)**0.5
                            'base_function': lambda x:  x, 'peaks_info': [[110, 2.95, 4000], [130, 3, 4000], [180, 3.41, 8000], [230, 3.75, 3500],
                                                                          [285, 4.14, 5000], [335, 4.4, 5000], [375, 4.7, 2000],
                                                                          [410, 4.95, 4800], [450, 5.2, 3200]]},
    'single_peaks_biggest': {'base_intensity': 100, 'base_amplitude': 100, 'label': 'synthesized_singlets', 'length': 500,  # with FWHM cal: 2+0.03*(E+0.05*E**2)**0.5
                            'base_function': lambda x:  x, 'peaks_info': [[110, 2.95, 12000], [130, 3, 8000], [180, 3.41, 6000], [230, 3.75, 7000],
                                                                          [285, 4.14, 10000], [335, 4.4, 10000], [375, 4.7, 5000],
                                                                          [410, 4.95, 7800], [450, 5.2, 6200]]},
    'test': {'base_intensity': 100, 'base_amplitude': 100, 'label': 'synthesized_singlets', 'length': 500,  # with FWHM cal: 2+0.03*(E+0.05*E**2)**0.5
                            'base_function': lambda x:  x, 'peaks_info': [[130, 3, 4000], [180, 3.41, 8000], [230, 3.75, 3500],
                                                                          [285, 4.14, 5000], [335, 4.4, 5000], [375, 4.7, 2000],
                                                                          [410, 4.95, 4800], [450, 5.2, 3200],
                                                                          [30, 2.33, 2000], [45, 2.42, 2000], [70, 2.57, 4000], [77, 2.61, 1500],
                                                                    [110, 2.83, 1500], [117, 2.87, 500], [126, 2.90, 2000],
                                                                    [150, 3.07, 2400], [165, 3.16, 1600],
                                                                    [275, 3.83, 2000], [290, 3.93, 2500], [305, 4.02, 1500],
                                                                    [310, 4.04, 1500], [317, 4.1, 500], [344, 4.22, 2000], 
                                                                    [470, 5.01, 2600]]},
}


for peak_info in simuspecs['test']['peaks_info']:
    peak_info[0] *= 2
    peak_info[2] *= 2
simuspecs['test']['length'] *= 2