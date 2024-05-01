from tokenize import blank_re
import numpy as np
import matplotlib.pyplot as plt
import xml.etree.ElementTree as ET
import xml.dom.minidom as minidom
import pickle
import pandas as pd

from time import strftime, localtime
from pathlib import Path
from copy import deepcopy
from typing import Callable, Any
from numpy.typing import NDArray
from xml.dom.minidom import parse as parse_xml
from re import match

from .plot_settings import colors
from .PeakRegion import Region, Calibration
from .globals import CALIBRATION_METHOD


class Spectrum(np.ndarray):
    """
    Spetrum is the msot basic object of GAMUT. It is used to store data of a gamma spectrum.
    It is a view of 1-d numpy.ndarray.
    
    Attributes:
        label: str, label of the spectrum.
        ergcal: Calibration, calibration of the energy axis.
        FWHMcal: Calibration, calibration of the FWHM axis.
        attrs: dict, other attributes.
        regions: list[Region], regions of interest.
    
    Methods:
        __init__: Initialize the Spectrum object.
        __array_finalize__: Finalize the Spectrum object.
        __getitem__: Get item from the Spectrum object.
        __rshift__: Right shift operator.
        length: Get the length of the Spectrum object.
        indexes: Get the indexes of the Spectrum object.
        energies: Get the energies of the Spectrum object.
        counts: Get the counts of the Spectrum object.
        npeaks: Get the number of peaks of the Spectrum object.
        peaks: Get the peaks of the Spectrum object.
        covariance_propagation_vector: Get the covariance propagation vector of the Spectrum object.
        propagate_covariance: Propagate the covariance of the Spectrum object.
        covariance: Get the covariance of the Spectrum object.
        slice: Slice the Spectrum object.
        copy: Copy the Spectrum object.
        estimate_area: Estimate the area of the Spectrum object.
        estimate_fitness: Estimate the fitness of the Spectrum object.
        plot: Plot the Spectrum object.
        plot_regions: Plot the regions of the Spectrum object.
        plot_peaks: Plot the peaks of the Spectrum object.
        export_to_xml: Export the Spectrum object to xml file.
        from_GammaVision: Import the Spectrum object from GammaVision file.
        from_pickle: Import the Spectrum object from pickle file.
        export_to_pickle: Export the Spectrum object to pickle file.
        export_to_GammaVision: Export the Spectrum object to GammaVision file.
        export_to_pandas: Export the Spectrum object to pandas DataFrame.
        from_excel: Import the Spectrum object from excel file.
        from_xml: Import the Spectrum object from xml file.
    """
    def __new__(cls, counts: list | NDArray[np.float64], *args, **kargs):
        self = np.asarray(counts, dtype=np.float64).view(cls)
        self._covariance_propagation_vector = np.array([1], dtype=np.float64)
        return self

    def __array_finalize__(self, obj):
        if obj is None:
            return
        self.label = deepcopy(getattr(obj, 'label', 'Spectrum'))
        self.ergcal = deepcopy(getattr(obj, 'ergcal', Calibration()))
        self.FWHMcal = deepcopy(getattr(obj, 'FWHMcal', Calibration()))
        self._covariance_propagation_vector = deepcopy(getattr(obj, '_covariance_propagation_vector', np.array([1], dtype=np.float64)))
        self.attrs = deepcopy(getattr(obj, 'attrs', {}))
        self.regions = deepcopy(getattr(obj, 'regions', []))

    def __init__(self, counts: list | np.ndarray, label: str | None = None,
                 regions: list[Region] | None = None, ergcal: Calibration | None = None,
                 FWHMcal: Calibration | None = None, attrs: dict | None = None):
        
        # input legimacy check
        if label is None:
            label = 'spectrum'
        if regions is None:
            regions = []
        if ergcal is None:
            ergcal = Calibration()  # default linear calibration
        if FWHMcal is None:
            FWHMcal = Calibration()
        if attrs is None:
            attrs = {}

        self.label = label
        self.ergcal = ergcal
        self.FWHMcal = FWHMcal
        self.attrs = attrs
        self.regions = regions

    def __array_wrap__(self, out_arr, context=None):
        if out_arr.ndim == 0:
            return out_arr.item()
        else:
            return np.ndarray.__array_wrap__(self, out_arr, context)

    def __repr__(self):
        return f"Spectrum[{self.label}]{self.shape}"

    def __str__(self):
        return f"Spectrum[{self.label}]"

    def __reduce__(self):
        pickled_state = super().__reduce__()
        if isinstance(pickled_state[2], tuple):
            new_state = pickled_state[2] + (self.label, self.regions, self.ergcal, self.FWHMcal, self.attrs)
            return (pickled_state[0], pickled_state[1], new_state)
        else:
            raise ValueError("pickled_state[2] is not a tuple(Unsolved problem)")

    def __setstate__(self, state):
        self.label = state[-5]
        self.regions = state[-4]
        self.ergcal = state[-3]
        self.FWHMcal = state[-2]
        self.attrs = state[-1]
        super().__setstate__(state[0:-5])

    def __getitem__(self, key):
        sliced = super().__getitem__(key)
        if not isinstance(key, slice):
            return sliced

        start = key.start if key.start is not None else 0
        stop = key.stop if key.stop is not None else self.length

        if self.ergcal.by == 'data':
            data = self.ergcal.data.copy()
            data[:, 0] -= start
            ergcal = Calibration(method=self.ergcal.method, data=data)
        else:
            params = self.ergcal.params.copy()
            params[0] += self.ergcal(start)
            ergcal = Calibration(method=self.ergcal.method, params=params)

        regions = []
        for region in self.regions:
            if (region.right <= start) or (region.left >= stop):
                continue
            else:
                region_copy = region.copy()
                region_copy.left = max(region_copy.left-start, 0)
                region_copy.right = min(region_copy.right-start, stop-start)
                for peak in region_copy.peaks:
                    peak['location'] = min(max(peak['location']-start, 0), stop-start)
                    if 'center' in peak.keys():
                        peak['center'] = min(max(peak['center']-start, 0), stop-start)
                regions.append(region_copy)

        return Spectrum(sliced, label=self.label, regions=regions, ergcal=ergcal, FWHMcal=self.FWHMcal, attrs=self.attrs)

    def __rshift__(self, operator: Callable):
        return operator(self)

    @property
    def length(self):
        return self.shape[0]

    @property
    def indexes(self):
        return np.arange(self.length)

    @property
    def energies(self) -> np.ndarray:
        return self.ergcal(self.indexes)

    @property
    def counts(self):
        return np.asarray(self)

    @property
    def npeaks(self):
        return sum([len(region.peaks) for region in self.regions])

    @property
    def peaks(self):
        return [peak for region in self.regions for peak in region.peaks]

    @property
    def covariance_propagation_vector(self) -> np.ndarray:
        return self._covariance_propagation_vector

    def propagate_covariance(self, covariance_propagation_vector: np.ndarray) -> None:
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

    def slice(self, erg: tuple) -> "Spectrum":
        ergstart = erg[0] if erg[0] is not None else self.energies.min()
        ergstop = erg[1] if erg[1] is not None else self.energies.max()
        start = np.abs(self.energies - ergstart).argmin()
        stop = np.abs(self.energies - ergstop).argmin()
        return self[start:stop+1].view(Spectrum)

    def copy(self):
        return deepcopy(self)

    def estimate_area(self, region: Region) -> tuple:
        """
        Estimate peak area of spectrum in a region based on classic counting method and linear baseline assumption.
        """
        total = sum(self[region.indexes])
        len_left = min(2, region.left)
        len_right = min(2, self.length - region.right - 1)
        baseline = (self[region.left-len_left: region.left+1].mean() + self[region.right: region.right+len_right+1].mean()) * region.length / 2
        return total, baseline

    def estimate_fitness(self, region: Region) -> float:
        """
        Estimate fitness of spectrum in a region based on fitted parameters.
        """
        fcounts = region.fit_peaks() + region.fit_baseline()
        return 1 - ((fcounts - self[region.indexes])**2).sum() / (((self[region.indexes] - self[region.indexes].mean())**2).sum() + 1)

    def plot(self, ax: plt.Axes | None = None, chinese_label: bool = False, *args, **kargs) -> plt.Axes:
        """
        Plot the spectrum.
        """
        if ax is None:
            ax: plt.Axes = plt.gca()
        ax.plot(self.energies, self, *args, **kargs, label=self.label)
        ax.set_ylim(0, max(self[10:])*1.5)
        ax.set_xlim(self.energies.min(), self.energies.max())
        if chinese_label:
            ax.set_xlabel('能量/keV')
            ax.set_ylabel('计数')
        else:
            ax.set_xlabel('Energy/keV')
            ax.set_ylabel('Counts')
        return ax

    def plot_regions(self, *args, ax: plt.Axes | None = None, chinese_label: bool = False, **kargs) -> plt.Axes:
        """
        Plot regions of the spectrum, with peaks marked as stars.
        """
        if not hasattr(self, 'regions'):
            raise ValueError(f"{self} does not have regions. Call Gamut.PeakSearcher first.")
        if ax is None:
            ax: plt.Axes = plt.gca()
        for i, region in enumerate(self.regions):
            ax.plot(self.ergcal(region.indexes), self[region.indexes], *args, **kargs, color=colors(i), linewidth=2, alpha=0.6)
            for peak in region.peaks:
                ax.plot(self.ergcal(peak['location']), self[peak['location']], *args, **kargs, color=colors(i), marker='v', markersize=6)
        if chinese_label:
            ax.set_xlabel('能量/keV')
            ax.set_ylabel('计数')
        else:
            ax.set_xlabel('Energy/keV')
            ax.set_ylabel('Counts')
        return ax

    @staticmethod
    def _gaussian(indexes, center, stderr):
        return np.exp(-(indexes - center) ** 2 / (2 * stderr ** 2))

    # def fit_spectrum(self, baseline, min_fitness: float | None = None) -> "Spectrum":
    #     if baseline is not None:
    #         baseline = baseline.copy().astype(np.float64)
    #     else:
    #         baseline = np.zeros(self.length)
    #     no_baseline = Spectrum(self - baseline)
    #     fitted_spectrum = baseline.copy()
    #     for region in self.regions:
    #         if (min_fitness is None) or (no_baseline.estimate_fitness(region) > min_fitness):
    #             fitted_baseline = region.fit_baseline()
    #             baseline[region.indexes] += fitted_baseline
    #             fitted_spectrum[region.indexes] += fitted_baseline
    #             for peak in region.peaks:
    #                 if ('stderr' in peak.keys()) and ('center' in peak.keys()):
    #                     fitted_indexes = region.indexes
    #                     fitted_peak = region.fit_peak(peak, fitted_indexes)
    #                     fitted_spectrum[fitted_indexes] += fitted_peak
    #     return fitted_spectrum

    def plot_peaks(self, min_fitness: float | None = None, axes: plt.Axes | None = None, chinese_label: bool = False,
                   baseline: np.ndarray | None = None, plot_baseline: bool = True, plot_fitted: bool = True, *args, **kargs) -> plt.Axes:
        """
        Plot detailed fitted peaks in regions of the spectrum.
        Baseline is optional to reconstruct the unstripped spectrum.
        """
        if baseline is not None:
            baseline = baseline.copy().astype(np.float64)
        else:
            baseline = np.zeros(self.length)
        no_baseline = Spectrum(self - baseline)

        fitted_spectrum = baseline.copy()
        if axes is None:
            axes: plt.Axes = plt.gca()

        for region in self.regions:
            if (min_fitness is None) or (no_baseline.estimate_fitness(region) > min_fitness):
                fitted_baseline = region.fit_baseline()
                baseline[region.indexes] += fitted_baseline
                fitted_spectrum[region.indexes] += fitted_baseline
                for peak in region.peaks:
                    if ('stderr' in peak.keys()) and ('center' in peak.keys()):
                        fitted_indexes = region.indexes
                        fitted_peak = region.fit_peak(peak, fitted_indexes)
                        axes.fill_between(x=self.ergcal(fitted_indexes),
                                          y1=baseline[fitted_indexes], y2=baseline[fitted_indexes] + fitted_peak, alpha=0.5)
                        fitted_spectrum[fitted_indexes] += fitted_peak
        if plot_baseline:
            axes.plot(self.energies, baseline, alpha=0.4, linewidth=1, label=self.label+'(baseline)', linestyle='-')
        if plot_fitted:
            axes.plot(self.energies, fitted_spectrum, linewidth=1, label=self.label+'(fitted)', linestyle='-')
        if chinese_label:
            axes.set_xlabel('能量/keV')
            axes.set_ylabel('计数')
        else:
            axes.set_xlabel('Energy/keV')
            axes.set_ylabel('Counts')
        return axes, fitted_spectrum

    def export_to_xml(self, filepath: str | Path,
                      data: dict[str, bool] = {'ergcal': True, 'ergcal.verbose': False,
                                               'FWHMcal': True, 'FWHMcal.verbose': False,
                                               'regions': True, 'regions.verbose': False,
                                               'regions.peaks': True, 'regions.peaks.verbose': False}):
        """
        Export the report to a xml file.
        """
        spec = ET.Element('Spectrum')
        spec.attrib['name'] = self.label
        spec.attrib['length'] = str(self.length)

        if data['ergcal']:
            ergcal_el = ET.SubElement(spec, 'ergcal')
            ergcal_el.attrib['method'] = self.ergcal.method
            ergcal_el.attrib['params'] = " ".join([f"{val:.4f}" for val in self.ergcal.params])
            if data['ergcal.verbose'] and hasattr(self.ergcal, 'data'):
                ET.SubElement(spec, 'data_indexes').text = " ".join([f"{val:.2f}" for val in self.ergcal.data_indexes])
                ET.SubElement(spec, 'data_values').text = " ".join([f"{val:.2f}" for val in self.ergcal.data_values])
                ergcal_el.attrib['fitness'] = f"{self.ergcal.fitness:.2f}"

        if data['FWHMcal']:
            FWHMcal_el = ET.SubElement(spec, 'FWHMcal')
            FWHMcal_el.attrib['method'] = self.FWHMcal.method
            FWHMcal_el.attrib['params'] = " ".join([f"{val:.4f}" for val in self.FWHMcal.params])
            if data['FWHMcal.verbose'] and hasattr(self.FWHMcal, 'data'):
                ET.SubElement(spec, 'data_ind').text = " ".join([f"{val:.2f}" for val in self.FWHMcal.data[0, :]])
                ET.SubElement(spec, 'data_val').text = " ".join([f"{val:.2f}" for val in self.FWHMcal.data[1, :]])
                FWHMcal_el.attrib['fitness'] = f"{self.FWHMcal.fitness:.2f}"

        if data['regions']:
            regions_el = ET.SubElement(spec, 'regions')
            for region in self.regions:
                region_el = ET.SubElement(regions_el, 'region')
                region_el.attrib['N_peaks'] = str(len(region.peaks))
                region_el.attrib['ind_left'] = str(region.left)
                region_el.attrib['ind_right'] = str(region.right)

                if data['regions.verbose']:
                    region_el.attrib['erg_left'] = f"{self.ergcal(region.left):.2f}"
                    region_el.attrib['erg_right'] = f"{self.ergcal(region.right):.2f}"
                    try:
                        region_el.attrib['fitness'] = f"{self.estimate_fitness(region):.2f}"
                        region_el.attrib['slope'] = f"{region.slope:.2f}"
                        region_el.attrib['offset'] = f"{region.offset:.2f}"
                    except:
                        pass
                
                if data['regions.peaks']:
                    for peak in region.peaks:
                        peak_el = ET.SubElement(region_el, 'peak')
                        peak_el.attrib['energy'] = f"{self.ergcal(peak['location']):.2f}"
                        for key, val in peak.items():
                            peak_el.attrib[key] = str(val) if isinstance(val, int) else f"{val:.2f}"

        ET.SubElement(spec, 'counts').text = " ".join([f"{val:.2f}" for val in self.counts])

        tree = ET.ElementTree(spec)
        tree.write(filepath, encoding='utf-8', xml_declaration=True)
        dom = parse_xml(str(filepath))
        pdom = dom.toprettyxml(indent='\t', newl='\n')
        dom.unlink()
        with open(filepath, 'w') as fileopen:
            fileopen.write(pdom)

        # spec = ET.Element('spectrum')
        # spec.attrib['name'] = self.label
        # spec.attrib['length'] = str(self.length)

        # if data['ergcal'] and hasattr(self, 'ergcal'):
        #     ergcal_el = ET.SubElement(spec, 'ergcal')
        #     ergcal_el.attrib['method'] = self.ergcal.method
        #     ergcal_el.attrib['params'] = " ".join([f"{val:.4f}" for val in self.ergcal.params])
        #     if hasattr(self.ergcal, 'data'):
        #         ET.SubElement(spec, 'data_indexes').text = " ".join([f"{val:.2f}" for val in self.ergcal.data[0, :]])
        #         ET.SubElement(spec, 'data_values').text = " ".join([f"{val:.2f}" for val in self.ergcal.data[1, :]])

        # if hasattr(self, 'FWHMcal'):
        #     FWHMcal_el = ET.SubElement(spec, 'FWHMcal')
        #     FWHMcal_el.attrib['method'] = self.FWHMcal.method
        #     FWHMcal_el.attrib['params'] = " ".join([f"{val:.4f}" for val in self.FWHMcal.params])
        #     if hasattr(self.FWHMcal, 'data'):
        #         ET.SubElement(spec, 'data_ind').text = " ".join([f"{val:.2f}" for val in self.FWHMcal.data[0, :]])
        #         ET.SubElement(spec, 'data_val').text = " ".join([f"{val:.2f}" for val in self.FWHMcal.data[1, :]])
        #     if hasattr(self.FWHMcal, '_fitness'):
        #         FWHMcal_el.attrib['fitness'] = f"{self.FWHMcal._fitness:.2f}"

        # if hasattr(self, 'regions'):
        #     regions_el = ET.SubElement(spec, 'regions')
        #     for region in self.regions:
        #         region_el = ET.SubElement(regions_el, 'region')
        #         region_el.attrib['ind_left'] = str(region.left)
        #         region_el.attrib['ind_right'] = str(region.right)
        #         region_el.attrib['erg_left'] = f"{self.ergcal(region.left):.2f}"
        #         region_el.attrib['erg_right'] = f"{self.ergcal(region.right):.2f}"
        #         region_el.attrib['N_peaks'] = str(len(region.peaks))
        #         try:
        #             region_el.attrib['fitness'] = f"{self.estimate_fitness(region):.2f}"
        #             region_el.attrib['slope'] = f"{region.slope:.2f}"
        #             region_el.attrib['offset'] = f"{region.offset:.2f}"
        #         except:
        #             pass
        #         for peak in region.peaks:
        #             peak_el = ET.SubElement(region_el, 'peak')
        #             for key, val in peak.items():
        #                 peak_el.attrib[key] = str(val) if isinstance(val, int) else f"{val:.2f}"
        #             peak_el.attrib['energy'] = f"{self.ergcal(peak['location']):.2f}"

        # ET.SubElement(spec, 'counts').text = " ".join([f"{val:.2f}" for val in self.counts])

        # tree = ET.ElementTree(spec)
        # tree.write(filepath, encoding='utf-8', xml_declaration=True)
        # dom = parse_xml(str(filepath))
        # pdom = dom.toprettyxml(indent='\t', newl='\n')
        # dom.unlink()
        # with open(filepath, 'w') as fileopen:
        #     fileopen.write(pdom)

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
        if data_index is None:
            raise ValueError(f"{filename} is not a valid GammaVision file.")
        index = data_index + 2
        counts = []
        while filelines[index].strip().isdigit():
            counts.append(int(filelines[index].strip()))
            index += 1
            if index >= len(filelines):
                break
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
            self.ergcal = Calibration(method='linear', params=[float(incerp), float(slope)])
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

    def export_to_pandas(self):
        """
        Export to excel file.
        """
        results = []
        for region in self.regions:
            for peak in region.peaks:

                val = {'left': region.left, 'right': region.right,
                        'energy': self.ergcal(peak['location'])}
                attrs = ['location', 'height', 'stderr', 'area', 'sig_area2', 'fitness']
                val.update(dict([(attr, peak[attr]) for attr in attrs if attr in peak.keys()]))
                if 'stderr' in peak.keys():
                    val.update({'fitness': self.estimate_fitness(region)})  

                results.append(val)
        results = pd.DataFrame(results)
        return results
        
    @classmethod
    def from_excel(cls, filename: str | Path, counts_column: str, energy_column: str | None = None):
        """
        Import spectrum from columns of excel files.
        """
        df = pd.read_excel(str(filename), index_col=0)
        if counts_column not in df.columns:
            raise ValueError(f"{counts_column} not in df.columns")
        else:
            self = cls(df[counts_column].to_numpy())
            if energy_column and energy_column not in df.columns:
                raise ValueError(f"{energy_column} not in df.columns")
            else:
                energies = df[energy_column].to_numpy()
                data = list(np.vstack((np.arange(energies.shape[0]), energies)))
                self.ergcal = Calibration(method='linear', data=data)
        return self

    @classmethod
    def from_xml(cls, filename: str | Path):
        """
        Import spectrum from xml files.
        """
        spec_xml = parse_xml(str(filename))
        if (spec_el := spec_xml.getElementsByTagName('counts')[0]) is None:
            raise ValueError(f"{filename} is not a valid xml file containing spectrum counts.")
        else:
            self = cls(counts=[float(val) for val in spec_el.firstChild.nodeValue.split()],
                    label=spec_el.getAttribute('name'))

        if (ergcal_el := spec_xml.getElementsByTagName('ergcal')[0]) is not None:
            self.ergcal = Calibration(method = ergcal_el.getAttribute('method'),
                                      params=[float(val) for val in ergcal_el.getAttribute('params').split()])

        if (FWHMcal_el := spec_xml.getElementsByTagName('FWHMcal')[0]) is not None:
            self.FWHMcal = Calibration(method=FWHMcal_el.getAttribute('method'),
                                       params=[float(val) for val in FWHMcal_el.getAttribute('params').split()])

        if (regions_el := spec_xml.getElementsByTagName('regions')[0]) is not None:
            for region_el in regions_el.getElementsByTagName('region'):
                region = Region(int(region_el.getAttribute('ind_left')), int(region_el.getAttribute('ind_right')))

                try:
                    region.slope = float(region_el.getAttribute('slope'))
                    region.offset = float(region_el.getAttribute('offset'))
                except:
                    pass

                for peak_el in region_el.getElementsByTagName('peak'):
                    peak = dict([(attr.name, float(attr.value)) for attr in peak_el.attributes.values()])
                    peak['location'] = int(peak['location'])
                    region.peaks.append(peak)
                    
                self.regions.append(region)

        return self

    @classmethod
    def from_MCNP(cls, filepath: str, tally_id: int = 8):
            
        with open(filepath, 'r') as fileopen:
            filelines = fileopen.readlines()   

        index = 0
        while f'tally type {tally_id}' not in filelines[index]:
            index += 1
        index += 1

        energies, counts = [], []
        while 'energy' not in filelines[index]:
            index += 1
        index += 1
        while 'total' not in filelines[index]:
            erg, eff, err = filelines[index].split()
            energies.append(float(erg))
            counts.append(float(eff))
            index += 1

        while 'run terminated when' not in filelines[index]:
            index += 1
        match_results = match(r'\s*run terminated when \s* ([0-9]*) \s* particle histories were done.', filelines[index])
        if match_results:
            NPS = int(match_results.group(1))
        else:
            NPS = 1

        counts = np.array(counts) * NPS
        energies = np.array(energies) * 1000

        return cls(counts, ergcal=Calibration(method='linear', data=np.stack([np.arange(len(counts)), energies], axis=1)))


class SimulatedSpectrum(Spectrum):
    """
    Simulated Spectrum Generator for algorithm test.
    """
    def __new__(cls, length: int = 100, *args, **kargs):
        counts = np.zeros(length)
        self = super().__new__(cls, counts)
        return self

    def __init__(self,
                 length: int = 100,
                 peaks_info: list = [[50, 3, 2000]],
                 pack_peak_info: bool = True,
                 base_intensity: int = 0,
                 base_amplitude: int = 0,
                 base_function: Callable = lambda x: x,
                 base_noised: bool = True,
                 random_seed: int = 0,
                 label: str | None = None,
                 ergcal: Calibration | None = None,
                 FWHMcal: Calibration | None = None
                 ):

        self.regions_info = peaks_info
        self.base_intensity = base_intensity
        self.base_amplitude = base_amplitude
        self.base_function = base_function
        self.base_noised = base_noised
        self.random_seed = random_seed
        np.random.seed(self.random_seed)

        self.baseline = self._gen_base()
        self.peakform = self._gen_peaks()
        self[:] = self.baseline[:] + self.peakform[:]
        if self.base_noised:
            self[:] = np.random.normal(loc=self[:], scale=self[:]**0.5, size=(self.length,))
        super().__init__(self.baseline + self.peakform, label, ergcal=ergcal, FWHMcal=FWHMcal)
        
        if pack_peak_info:
            self.regions = []
            for peak_info in peaks_info:
                left, right = int(peak_info[0]-peak_info[1]*3), int(peak_info[0]+peak_info[1]*3)
                region = Region(left, right,
                                peaks=[{'location': round(peak_info[0]), 
                                        'center': peak_info[0], 
                                        'stderr': peak_info[1], 
                                        'area': peak_info[2],
                                        'height': peak_info[2] / peak_info[1] / (2*np.pi)**0.5}])
                region.slope = 0
                region.offset = self[round(region.peaks[0]['center'])] - region.peaks[0]['height']
                self.regions.append(region)

    def _gen_base(self):
        base = self.base_function(np.linspace(0, 1, self.length))
        base = (base - base.min()) / (base.ptp() + 1E-3)
        base = base * self.base_amplitude + self.base_intensity
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