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

from ..globals import colors
from .Spectrum import Spectrum


class SpectrumPlotter:
    
    def plot(self, ax: plt.Axes | None = None, chinese_label: bool = False, *args, **kargs) -> plt.Axes:
        """
        Plot the spectrum.

        :param ax: The axes to plot on.
        :param chinese_label: Whether to use Chinese labels.
        :param args: Other arguments to pass to plt.plot
        :param kargs: Other keyword arguments to pass to plt.plot
        
        :return: The axes object.
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

    def plot_regions(self, ax: plt.Axes | None = None, chinese_label: bool = False, *args, **kargs) -> plt.Axes:
        """
        Plot regions of the spectrum, with peaks marked as stars.
        
        :param ax: The axes to plot on.
        :param chinese_label: Whether to use Chinese labels.
        :param args: Other arguments to pass to plt.plot
        :param kargs: Other keyword arguments to pass to plt.plot
        
        :return: The axes object.
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
        
        :param min_fitness: The minimum fitness of the region to be fitted.
        :param axes: The axes to plot on.
        :param chinese_label: Whether to use Chinese labels.
        :param baseline: The baseline to be fitted.
        :param plot_baseline: Whether to plot the baseline.
        :param plot_fitted: Whether to plot the fitted spectrum.
        :param args: Other arguments to pass to plt.plot
        :param kargs: Other keyword arguments to pass to plt.plot
        
        :return: The axes object.
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