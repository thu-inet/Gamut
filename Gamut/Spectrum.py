import numpy as np
import matplotlib.pyplot as plt
import xml.etree.ElementTree as ET
import pickle
import pandas as pd

from time import strftime, localtime
from pathlib import Path
from copy import deepcopy
from typing import Callable
from xml.dom.minidom import parse as parse_xml

from .plot_settings import colors
from .PeakRegion import Region, Calibration


class Spectrum(np.ndarray):
    """
    1. The basic object of Gamut, the Spectrum. It is used to store data of a spectrum.
    2. Intrinsiclly it is a view of 1-d numpy.ndarray, so its all methods and attributes are available.
    3. Contains three main sub-objects: regions, ergcal, and FWHMcal. Attrs to store other data.
    """
    def __new__(cls, counts: np.ndarray, *args, **kargs):
        self = np.asarray(counts, dtype=np.float32).view(cls)
        self._covariance_propagation_vector = np.array([1], dtype=np.float32)
        return self

    def __array_finalize__(self, obj):
        if obj is None:
            return
        self._label = deepcopy(getattr(obj, '_label', 'Spectrum'))
        self._regions = deepcopy(getattr(obj, '_regions', []))
        self._ergcal = deepcopy(getattr(obj, '_ergcal', Calibration()))
        self._FWHMcal = deepcopy(getattr(obj, '_FWHMcal', Calibration()))
        self._covariance_propagation_vector = deepcopy(getattr(obj, '_covariance_propagation_vector', np.array([1], dtype=np.float32)))
        self._attrs = deepcopy(getattr(obj, '_attrs', None))

    def __init__(self, counts: list | np.ndarray, label: str = None,
                 regions: list[Region] = None, ergcal: Calibration = None,
                 FWHMcal: Calibration = None, attrs: dict = None):
        if label is None:
            label = 'Spectrum'
        if regions is None:
            regions = []
        if ergcal is None:
            ergcal = Calibration()
        if FWHMcal is None:
            FWHMcal = Calibration()
        if attrs is None:
            attrs = {}

        self._label = label
        self._regions = regions
        self._ergcal = ergcal
        self._FWHMcal = FWHMcal
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

    def __reduce__(self):
        pickled_state = super().__reduce__()
        new_state = pickled_state[2] + (self._label, self._regions, self._ergcal, self._FWHMcal, self._attrs)
        return (pickled_state[0], pickled_state[1], new_state)

    def __setstate__(self, state):
        self._label = state[-5]
        self._regions = state[-4]
        self._ergcal = state[-3]
        self._FWHMcal = state[-2]
        self._attrs = state[-1]
        super().__setstate__(state[0:-5])

    @property
    def label(self):
        return self._label

    @label.setter
    def label(self, value):
        self._label = value

    @property
    def energies(self):
        return self.ergcal(np.arange(len(self)))

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
    def regions(self):
        return self._regions

    @regions.setter
    def regions(self, regions: list[Region]):
        self._regions = regions

    @property
    def ergcal(self):
        return self._ergcal

    @ergcal.setter
    def ergcal(self, ergcal):
        self._ergcal = ergcal

    @property
    def FWHMcal(self):
        return self._FWHMcal

    @FWHMcal.setter
    def FWHMcal(self, FWHMcal):
        self._FWHMcal = FWHMcal

    @property
    def covariance_propagation_vector(self):
        return self._covariance_propagation_vector

    def propagate_covariance(self, covariance_propagation_vector):
        self._covariance_propagation_vector = np.convolve(covariance_propagation_vector, self.covariance_propagation_vector)

    def covariance(self, region: Region):
        propagation_hwidth = (self._covariance_propagation_vector.shape[0] - 1) // 2  # construct the propagation matrix using stride trick
        zero_pad = np.zeros(region.length - 1, dtype=self.covariance_propagation_vector.dtype)
        row_templet = np.concatenate((zero_pad, self.covariance_propagation_vector, zero_pad))
        stride = row_templet.strides[0]
        covariance_propagation_matrix = np.lib.stride_tricks.as_strided(row_templet[region.length - 1:], shape=(region.length, region.length + 2*propagation_hwidth), strides=(-stride, stride))
        left_width, right_width = min(region.left, propagation_hwidth), min(self.length - region.right - 1, propagation_hwidth)
        unpropagated_covariance = np.diag(self[region.left - left_width: region.right + right_width + 1])  # construct the unpropagated covariance matrix
        unpropagated_covariance = np.pad(unpropagated_covariance, (propagation_hwidth - left_width, propagation_hwidth - right_width), mode='edge')
        return np.einsum('ik,kg,jg->ij', covariance_propagation_matrix, unpropagated_covariance, covariance_propagation_matrix)

    def copy(self):
        return deepcopy(self)

    def area_estimate(self, region: Region) -> tuple:
        """
        Estimate peak area of spectrum in a region based on classic counting method and linear baseline assumption.
        """
        total = sum(self[region.indexes])
        baseline = (self[region.left] + self[region.right]) * region.length / 2
        return total, baseline

    def fitness_estimate(self, region: Region) -> float:
        """
        Estimate fitness of spectrum in a region based on fitted parameters.
        """
        fcounts = region.fit_peaks() + region.fit_baseline()
        return 1 - ((fcounts - self[region.indexes])**2).sum() / ((self[region.indexes] - self[region.indexes].mean())**2).sum()

    def plot(self, *args, axes: plt.Axes = None, **kargs) -> plt.Axes:
        """
        Plot the spectrum.
        """
        if axes is None:
            axes = plt.gca()
        axes.plot(self.energies, self, *args, **kargs, label=self._label)
        axes.set_ylim(0, )
        axes.set_xlim(0, )
        return axes

    def plot_regions(self, *args, axes: plt.Axes = None, **kargs) -> plt.Axes:
        """
        Plot regions of the spectrum, with peaks marked as stars.
        """
        if not hasattr(self, 'regions'):
            raise ValueError(f"{self} does not have regions. Call Gamut.PeakSearcher first.")
        if axes is None:
            axes = plt.gca()
        for i, region in enumerate(self.regions):
            axes.plot(self.ergcal(region.indexes), self[region.indexes], *args, **kargs, color=colors(i), linewidth=1, alpha=0.8)
            for peak in region.peaks:
                axes.plot(self.ergcal(peak['location']), self[peak['location']], *args, **kargs, color=colors(i), marker='*', markersize=8)
        return axes

    @staticmethod
    def _gaussian(indexes, center, stderr):
        return np.exp(-(indexes - center) ** 2 / (2 * stderr ** 2))

    def plot_peaks(self, *args, axes: plt.Axes = None, baseline: np.ndarray = None, **kargs) -> plt.Axes:
        """
        Plot detailed fitted peaks in regions of the spectrum.
        Baseline is optional to reconstruct the unstripped spectrum.
        """
        if baseline is not None:
            baseline = baseline.copy().astype(np.float64)
        else:
            baseline = np.zeros(self.length)
        fitted_spectrum = baseline.copy()
        if axes is None:
            axes = plt.gca()
        for region in self.regions:
            if self.fitness_estimate(region) > 0.1:
                fitted_baseline = region.fit_baseline()
                baseline[region.indexes] += fitted_baseline
                fitted_spectrum[region.indexes] += fitted_baseline
                for peak in region.peaks:
                    if 'stderr' in peak.keys():
                        # fitted_indexes = np.arange(max(int(peak['center']-4*peak['stderr']), 0), min(int(peak['center']+4*peak['stderr']), self.length-1))
                        fitted_indexes = region.indexes
                        fitted_peak = region.fit_peak(peak, fitted_indexes)
                        axes.fill_between(self.ergcal(fitted_indexes), baseline[fitted_indexes], baseline[fitted_indexes] + fitted_peak, alpha=0.5)
                        fitted_spectrum[fitted_indexes] += fitted_peak
        axes.plot(self.energies, baseline, alpha=0.4, linewidth=1, label=self.label+'(baseline)', linestyle='--')
        axes.plot(self.energies, fitted_spectrum, linewidth=1, label=self.label+'(fitted)')
        return axes

    def export_to_xml(self, filepath: str | Path):
        """
        Export the report to a xml file.
        """
        spec = ET.Element('spectrum')
        spec.attrib['name'] = self.label
        spec.attrib['length'] = str(self.length)

        if hasattr(self, 'ergcal'):
            ergcal_el = ET.SubElement(spec, 'ergcal')
            ergcal_el.attrib['method'] = self.ergcal.method
            ergcal_el.attrib['params'] = " ".join([f"{val:.4f}" for val in self.ergcal.params])
            if hasattr(self.ergcal, 'data'):
                ET.SubElement(spec, 'data_ind').text = " ".join([f"{val:.2f}" for val in self.ergcal.data[0, :]])
                ET.SubElement(spec, 'data_val').text = " ".join([f"{val:.2f}" for val in self.ergcal.data[1, :]])

        if hasattr(self, 'FWHMcal'):
            FWHMcal_el = ET.SubElement(spec, 'FWHMcal')
            FWHMcal_el.attrib['method'] = self.FWHMcal.method
            FWHMcal_el.attrib['params'] = " ".join([f"{val:.4f}" for val in self.FWHMcal.params])
            if hasattr(self.FWHMcal, 'data'):
                ET.SubElement(spec, 'data_ind').text = " ".join([f"{val:.2f}" for val in self.FWHMcal.data[0, :]])
                ET.SubElement(spec, 'data_val').text = " ".join([f"{val:.2f}" for val in self.FWHMcal.data[1, :]])
            if hasattr(self.FWHMcal, '_fitness'):
                FWHMcal_el.attrib['fitness'] = f"{self.FWHMcal._fitness:.2f}"

        if hasattr(self, 'regions'):
            regions_el = ET.SubElement(spec, 'regions')
            for region in self.regions:
                region_el = ET.SubElement(regions_el, 'region')
                region_el.attrib['ind_left'] = str(region.left)
                region_el.attrib['ind_right'] = str(region.right)
                region_el.attrib['erg_left'] = f"{self.ergcal(region.left):.2f}"
                region_el.attrib['erg_right'] = f"{self.ergcal(region.right):.2f}"
                region_el.attrib['N_peaks'] = str(len(region.peaks))
                region_el.attrib['fitness'] = f"{self.fitness_estimate(region):.2f}"
                region_el.attrib['slope'] = f"{region.slope:.2f}"
                region_el.attrib['offset'] = f"{region.offset:.2f}"
                for peak in region.peaks:
                    peak_el = ET.SubElement(region_el, 'peak')
                    for key, val in peak.items():
                        peak_el.attrib[key] = str(val) if isinstance(val, int) else f"{val:.2f}"
                    peak_el.attrib['energy'] = f"{self.ergcal(peak['location']):.2f}"

        ET.SubElement(spec, 'counts').text = " ".join([f"{val:.2f}" for val in self.counts])

        tree = ET.ElementTree(spec)
        tree.write(filepath, encoding='utf-8', xml_declaration=True)
        dom = parse_xml(filepath)
        pdom = dom.toprettyxml(indent='\t', newl='\n')
        dom.unlink()
        with open(filepath, 'w') as fileopen:
            fileopen.write(pdom)

    @classmethod
    def from_GammaVision(cls, filename: str | Path):
        """
        Import spectrum from Ortec .chn or .spe ASCII plain text file.
        """
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

        region_index = next((i for i, line in enumerate(filelines) if line.startswith('$ROI')), None)
        if region_index is not None:
            index = region_index + 2
            self.regions = []
            while filelines[index].strip().isdigit():
                left, right = filelines[index].split()
                self.regions.append(Region(int(left), int(right)))
            index += 1

        ergcal_index = next((i for i, line in enumerate(filelines) if line.startswith('$ENER_FIT')), None)
        if ergcal_index is not None:
            incerp, slope = filelines[ergcal_index+1].split()
            self.ergcal = Calibration(method='linear', params=[float(incerp), float(slope), 0])
        return self

    @classmethod
    def from_pickle(cls, filename: str | Path):
        """
        Import spectrum from pickled Spectrum file.
        """
        with open(filename, 'rb') as fopen:
            self = pickle.load(fopen)
        return self

    def export_to_pickle(self, filename: str | Path):
        """
        Export to binary pickle file.
        """
        with open(filename, 'wb') as fopen:
            pickle.dump(self, fopen)

    def export_to_GammaVision(self, filename: str | Path):
        """
        Export to Ortec plain text file.
        """
        with open(filename, 'w') as fopen:
            fopen.write("$SPEC_ID:\n\n")
            fopen.write("$SPEC_REM:\n\n")
            fopen.write("$DATE_MEA:\n" + strftime(r"%m/%d/%Y %H:%M:%S", localtime()) + "\n")
            fopen.write("$MEAS_TIME:\n" + "1000 1000\n")
            fopen.write(f"$DATA:\n0 {self.length-1}\n")
            fopen.write("\n".join([f"{round(c):>8d}" for c in self.counts]) + "\n")
            if len(self.regions) != 0:
                fopen.write(f"$ROI:\n{len(self.regions)}\n" + "\n".join([f"{region.left} {region.right}" for region in self.regions]))
            fopen.write("$PRESETS:\nNone\n0\n0\n")
            fopen.write(f"$ENER_FIT:\n{self.ergcal.params[0]:8.6f} {self.ergcal.params[1]:8.6f}\n")

    @classmethod
    def from_excel(cls, filename: str | Path, energy_column: str, counts_column: str):
        """
        Import spectrum from columns of excel files.
        """
        df = pd.read_excel(filename, index_col=0)
        self = cls(df[counts_column])
        self.ergcal = Calibration(method='linear', data=[np.arange(len(df[energy_column])), df[energy_column]])
        return self

    @classmethod
    def from_xml(cls, filename: str | Path):
        """
        Import spectrum from xml files.
        """
        spec_xml = parse_xml(str(filename))
        spec_el = spec_xml.getElementsByTagName('counts')[0]
        self = cls(counts=[float(val) for val in spec_el.firstChild.nodeValue.split()],
                   label=spec_el.getAttribute('name'))

        if (ergcal_el := spec_xml.getElementsByTagName('ergcal')[0]) is not None:
            self.ergcal = Calibration(method=ergcal_el.getAttribute('method'),
                                      params=[float(val) for val in ergcal_el.getAttribute('params').split()])

        if (FWHMcal_el := spec_xml.getElementsByTagName('FWHMcal')[0]) is not None:
            self.FWHMcal = Calibration(method=FWHMcal_el.getAttribute('method'),
                                       params=[float(val) for val in FWHMcal_el.getAttribute('params').split()])

        if (regions_el := spec_xml.getElementsByTagName('regions')[0]) is not None:
            for region_el in regions_el.getElementsByTagName('region'):
                region = Region(int(region_el.getAttribute('ind_left')), int(region_el.getAttribute('ind_right')))
                region.slope = float(region_el.getAttribute('slope'))
                region.offset = float(region_el.getAttribute('offset'))
                for peak_el in region_el.getElementsByTagName('peak'):
                    peak = dict([(attr.name, float(attr.value)) for attr in peak_el.attributes.values()])
                    peak['location'] = int(peak['location'])
                    region.peaks.append(peak)
                self.regions.append(region)
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
                 length: int = 100,
                 peaks_info: list = [[50, 3, 2000]],
                 base_intensity: int = 0,
                 base_amplitude: int = 0,
                 base_function: Callable = lambda x: x,
                 random_seed: int = 0,
                 label: str = None,
                 ergcal: Calibration = None,
                 FWHMcal: Calibration = None
                 ):

        self.regions_info = peaks_info
        self.base_intensity = base_intensity
        self.base_amplitude = base_amplitude
        self.base_function = base_function
        self.random_seed = random_seed
        np.random.seed(self.random_seed)

        self.predefined_regions = []
        for peak_info in peaks_info:
            self.predefined_regions.append(Region(int(peak_info[0]-peak_info[1]*3), int(peak_info[0]+peak_info[1]*3)))

        self.baseline = self._gen_base()
        self.peakform = self._gen_peaks()
        self[:] = self.baseline[:] + self.peakform[:]
        super().__init__(self.baseline + self.peakform, label, ergcal=ergcal, FWHMcal=FWHMcal)

    def _gen_base(self):
        base = self.base_function(np.linspace(0, 1, self.length))
        base = (base - base.min()) / (base.ptp() + 1E-3)
        base = base * self.base_amplitude + self.base_intensity
        base = np.random.normal(loc=base, scale=base**0.5, size=(self.length,))
        return base

    def _gen_peaks(self):
        peak = np.zeros((self.length,))
        for peak_info in self.regions_info:
            peak += self._gen_peak(*peak_info)
        return peak

    def _gen_peak(self, centroid, stderror, area):
        channels = np.arange(0, self.length)
        amplitude = area / stderror / (2*np.pi)**0.5
        # print(amplitude, stderror)
        return amplitude * np.exp(- (channels-centroid)**2 / stderror**2 / 2)


simuspecs = {
    'simple': {'base_intensity': 100, 'base_amplitude': 100, 'label': 'simple',  # FWHM~2.355*5~14.13, params[linear]=[0.14873, 0]
               'base_function': lambda x:  x**0 + 1, 'peaks_info': [[95, 6, 4000]], 'FWHMcal': Calibration('linear', params=[0, 0.14873, 0])},
    'double_slight': {'base_intensity': 100, 'base_amplitude': 100, 'label': 'double_slight',
                      'base_function': lambda x:  x, 'peaks_info': [[70, 6, 4000], [105, 7, 3000]], 'FWHMcal': Calibration('linear', params=[9.4225, 0.067285, 0])},
    'double_normal_narrow': {'base_intensity': 100, 'base_amplitude': 1000, 'label': 'double_normal_narrow',
                             'base_function': lambda x:  1 / (x + 0.1), 'peaks_info': [[90, 3, 4000], [100, 2.5, 3000]]},
    'double_normal': {'base_intensity': 100, 'base_amplitude': 100, 'label': 'double_normal',
                      'base_function': lambda x:  1 / (x + 1), 'peaks_info': [[70, 6, 4000], [100, 10, 3000]]},
    'double_severe': {'base_intensity': 200, 'base_amplitude': 100, 'label': 'double_severe',
                      'base_function': lambda x:  np.abs((x-0.4)), 'peaks_info': [[70, 7, 4000], [90, 10, 3000]]},
    'synthesized': {'base_intensity': 100, 'base_amplitude': 100, 'label': 'synthesized_multiplets',
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
                                                                          [410, 4.95, 7800], [450, 5.2, 6200]]}
}


