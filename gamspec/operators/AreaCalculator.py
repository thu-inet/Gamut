import numpy as np
from scipy.optimize import fmin, minimize, Bounds
from copy import deepcopy
from time import strftime, ctime
from typing import Literal, Callable

from ..Spectrum import Spectrum
from ..Operator import Operator
from ..PeakRegion import Region


class AverageAreaCalculator(Operator):

    def __init__(self, times: int,
                 label: str | None = None):
        """
        The most basic PeakAreaCalculator.

        :param base_mode: mode to strip counts of baseline.
        :param label: label for the operator
        """
        self._times = times
        if label is None:
            label = f"AreaCalculator"
        super().__init__(1, label)

    def __run__(self, spectra: list[Spectrum], *args, **kargs) -> Spectrum:
        calculated = spectra[0].copy()
        for region in calculated.regions:
            self._calculate(region, calculated)
        return calculated

    def _split_regions(self, region: Region, spectrum: Spectrum) -> list[int]:
        splited_region_indexes = [(spectrum[region.peaks[i]['location']: region.peaks[i+1]['location']]).argmin()+region.peaks[i]['location'] for i in range(region.npeaks-1)]
        splited_region_indexes = [region.left] + splited_region_indexes + [region.right+1]
        return splited_region_indexes

    def _calculate(self, region: Region, spectrum: Spectrum):
        if region.npeaks == 0:
            pass
        elif region.npeaks == 1:
            peak_areas = []
            for i in range(self._times):
                left = max(region.left - i, 0)
                right = min(region.right + i, len(spectrum)-1)
                total = spectrum[left: right+1].sum()
                baseline = (spectrum[left] + spectrum[right]) * (right-left+1) / 2
                peak_areas.append(total - baseline)
            region.peaks[0]['area'] = sum(peak_areas) / len(peak_areas)
        else:
            splited_region_indexes = self._split_regions(region, spectrum)
            for i, peak in enumerate(region.peaks):
                total = spectrum[splited_region_indexes[i]: splited_region_indexes[i+1]+1].sum()
                baseline = (spectrum[splited_region_indexes[i]] + spectrum[min(splited_region_indexes[i+1], len(spectrum)-1)]) * (splited_region_indexes[i+1] - splited_region_indexes[i] + 1) / 2
                peak['area'] = total - baseline


class AverageCovellPeakAreaCalculator(Operator):

    def __init__(self, init_width: int, times: int, label: str | None = None):
        self._init_width = init_width
        self._times = times
        if label is None:
            label = f"CovellCalculator"
        super().__init__(1, label)

    def __run__(self, spectra: list[Spectrum], *args, **kargs) -> Spectrum:
        calculated = spectra[0].copy()
        for region in calculated.regions:
            self._calculate(region, calculated)
        return calculated

    def _split_regions(self, region: Region, spectrum: Spectrum) -> list[int]:
        splited_region_indexes = [(spectrum[region.peaks[i]['location']: region.peaks[i+1]['location']]).argmin()+region.peaks[i]['location'] for i in range(region.npeaks-1)]
        splited_region_indexes = [region.left] + splited_region_indexes + [region.right+1]
        return splited_region_indexes

    def _calculate(self, region: Region, spectrum: Spectrum):
        if region.npeaks == 0:
            pass
        elif region.npeaks == 1:
            peak_areas = []
            for i in range(self._times):
                left = max(region.left - i, 0)
                right = min(region.right + i, len(spectrum)-1)
                total = spectrum[left: right+1].sum()
                baseline = (spectrum[left] + spectrum[right]) * (right-left+1) / 2
                peak_areas.append(total - baseline)
            region.peaks[0]['area'] = sum(peak_areas) / len(peak_areas)
        else:
            splited_region_indexes = self._split_regions(region, spectrum)
            for i, peak in enumerate(region.peaks):
                total = spectrum[splited_region_indexes[i]: splited_region_indexes[i+1]+1].sum()
                baseline = (spectrum[splited_region_indexes[i]] + spectrum[min(splited_region_indexes[i+1], len(spectrum)-1)]) * (splited_region_indexes[i+1] - splited_region_indexes[i] + 1) / 2
                peak['area'] = total - baseline


class RegionPeakFitter(Operator):

    def __init__(self, trial: int = 1, equal_width: bool = True,
                 maxiter: int = 100,
                 baseline: bool = True, variable_npeaks: bool = False, label: str | None = None):
        self._trial = trial
        self._maxiter = maxiter
        self._baseline = baseline
        self._equal_width = equal_width
        self._variable_npeaks = variable_npeaks
        if label is None:
            label = "RegionPeakFitter"
        super().__init__(1, label)

    def __run__(self, spectra: list[Spectrum], *args, **kwargs) -> Spectrum:
        fitted = spectra[0].copy()
        regions = [region for region in fitted.regions if region.npeaks > 0]
        for region in regions:
            print(f"Fitting Region: {region.left}~{region.right}, NPeaks={region.npeaks}, time={strftime(ctime())}")
            best_params, best_heights, _, _, npeaks, covariance_matrix = self._fit(region, fitted)
            if self._equal_width:
                region.peaks = []
                for i in range(npeaks):
                    if round(best_params[i]) in region.indexes:
                        peak = {}
                        peak['location'] = round(best_params[i])
                        peak['center'] = best_params[i]
                        peak['stderr'] = abs(best_params[-2])
                        peak['height'] = best_heights[i]
                        peak['area'] = peak['height'] * peak['stderr'] * (2*np.pi)**0.5

                        peak['sig_center'] = covariance_matrix[i, i]**0.5
                        peak['sig_stderr'] = covariance_matrix[npeaks, npeaks]**0.5
                        peak['sig_height'] = covariance_matrix[npeaks+2+i, npeaks+2+i]**0.5
                        expectance_stderr_height = covariance_matrix[npeaks+2+i, npeaks] + peak['stderr'] * peak['height']
                        expectance_stderr2_height2 = (peak['sig_stderr']**2 + peak['stderr']**2) * (peak['sig_height']**2 + peak['height']**2)
                        peak['sig_area'] = (2 * np.pi * (expectance_stderr2_height2 - expectance_stderr_height**2)) ** 0.5  - covariance_matrix[npeaks, npeaks+2+i] ** 2

                        cov_height_stderr = covariance_matrix[npeaks, npeaks+2+i]
                        peak['sig_area2'] = (2 * np.pi * (peak['height']**2 * peak['sig_stderr']**2 + peak['sig_height']**2 * peak['stderr']**2 + 2 * peak['height'] * peak['stderr'] * cov_height_stderr)) ** 0.5
                        region.peaks.append(peak)
            else:
                region.peaks = []
                for i in range(npeaks):
                    if round(best_params[2*i]) in region.indexes:
                        peak = {}
                        peak['location'] = round(best_params[2*i])
                        peak['center'] = best_params[2*i]
                        peak['stderr'] = abs(best_params[2*i+1])
                        peak['height'] = best_heights[i]
                        peak['area'] = peak['height'] * peak['stderr'] * (2*np.pi)**0.5

                        peak['sig_center'] = covariance_matrix[2*i, 2*i]**0.5
                        peak['sig_stderr'] = covariance_matrix[i*2+1, i*2+1]**0.5
                        peak['sig_height'] = covariance_matrix[2*npeaks+i+1, 2*npeaks+i+1]**0.5
                        expectance_stderr_height = covariance_matrix[2*npeaks+i+1, i*2+1] + peak['stderr'] * peak['height']
                        expectance_stderr2_height2 = (peak['sig_stderr']**2 + peak['stderr']**2) * (peak['sig_height']**2 + peak['height']**2)
                        peak['sig_area'] = (2 * np.pi * (expectance_stderr2_height2 - expectance_stderr_height**2)) ** 0.5 - covariance_matrix[i*2+1, 2*npeaks+i+1] ** 2

                        cov_height_stderr = covariance_matrix[i*2+1, 2*npeaks+i+1]
                        peak['sig_area2'] = (2 * np.pi * (peak['height']**2 * peak['sig_stderr']**2 + peak['sig_height']**2 * peak['stderr']**2 + 2 * peak['height'] * peak['stderr'] * cov_height_stderr)) ** 0.5
                        region.peaks.append(peak)

            print(f"Finish Fitting Region: time={strftime(ctime())}")
            region.slope = best_params[-1] * best_heights[-1]
            region.offset = best_heights[-1]
            region.covariance = covariance_matrix
            region.peaks = sorted(region.peaks, key=lambda peak: peak['location'])
        return fitted

    def _fit(self, region: Region, fitted: Spectrum):
        best_error = 1E8
        npeaks = region.npeaks
        heights = np.ones(npeaks + 1)  # mutable object to transfer data without return
        shapes = np.zeros((region.length, npeaks + 1))
        fcounts = fitted[region.indexes].copy()

        if self._equal_width:
            guesses = np.array([peak['location'] for peak in region.peaks])
            guesses = np.hstack((guesses, [region.length/npeaks/4, 0]))
        else:
            guesses = np.array([[peak['location'], region.length/npeaks/4] for peak in region.peaks]).ravel()
            guesses = np.hstack((guesses, 0))

        best_params, best_heights, best_shapes, best_fcounts = guesses.copy(), heights.copy(), shapes.copy(), fcounts.copy()
        for _ in range(self._trial):
            params, error = fmin(self._fitgaussian, x0=guesses,
                                 args=(region.indexes, fitted[region.indexes], shapes, heights, fcounts, npeaks),
                                 maxiter=self._maxiter,
                                 disp=False, full_output=True)[:2]

            if abs(error/best_error - 1) < 1E-2:
                break

            if (error <= best_error and heights[:-1].min() >= 0):
                best_error = error
                guesses = params * np.random.normal(1, 0.01, len(params))
                best_params = params.copy()
                best_heights = heights.copy()
                best_shapes = shapes.copy()
                best_fcounts = fcounts.copy()

        spectrum_covariance_matrix = fitted.covariance(region)
        covariance_matrix = self._covariance_estimate(region, npeaks, best_params, best_heights, spectrum_covariance_matrix)

        return best_params, best_heights, best_shapes, best_fcounts, npeaks, covariance_matrix

    def _fitgaussian(self, params: np.ndarray, indexes: np.ndarray, counts: np.ndarray,
                     shapes: np.ndarray, heights: np.ndarray, fcounts: np.ndarray, npeaks: int):
        if self._equal_width:
            for i in range(npeaks):
                shapes[:, i] = __class__._gaussian(indexes, params[i], params[-2])
        else:
            for i in range(npeaks):
                shapes[:, i] = __class__._gaussian(indexes, params[2*i], params[2*i+1])
        shapes[:, -1] = __class__._linear(indexes, params[-1]) * self._baseline
        heights[:] = np.linalg.lstsq(shapes, counts, rcond=None)[0]
        heights[:-1] = np.abs(heights[:-1])
        fcounts[:] = shapes @ heights
        return np.linalg.norm(fcounts - counts, ord=2)

    @staticmethod
    def _gaussian(x: np.ndarray, center: float, std: float):
        return np.exp(-(x - center) ** 2 / (2 * std ** 2))

    @staticmethod
    def _linear(x: np.ndarray, slope: float):
        return slope * x + 1

    def _gaussian_partial(self, indexes: np.ndarray, center: float, stderr: float, height: float):
        partials = np.zeros((len(indexes), 3))
        shape = np.exp(-(indexes - center) ** 2 / (2 * stderr**2))
        partials[:, 0] = height * shape * (indexes-center) / stderr**2
        partials[:, 1] = height * shape * (indexes - center) ** 2 / stderr**3
        partials[:, 2] = shape
        return partials

    def _covariance_estimate(self, region: Region, npeaks: int, best_params: np.ndarray, best_heights: np.ndarray, spectrum_covariance_matrix: np.ndarray) -> np.ndarray:
        partial_matrix = np.zeros((region.length, len(best_params)+len(best_heights)))
        if self._equal_width:
            for i in range(npeaks):
                partial = self._gaussian_partial(region.indexes, best_params[i], best_params[npeaks], best_heights[i])
                partial_matrix[:, i] = partial[:, 0].copy()
                partial_matrix[:, npeaks] += partial[:, 1].copy()
                partial_matrix[:, npeaks+2+i] = partial[:, 2].copy()
            partial_matrix[:, npeaks+1] = best_heights[-1] * region.indexes if self._baseline else 1
            partial_matrix[:, -1] = 1 + best_params[-1] * region.indexes
        else:
            for i in range(npeaks):
                partial = self._gaussian_partial(region.indexes, best_params[i*2], best_params[i*2+1], best_heights[i])
                partial_matrix[:, i*2: i*2+2] = partial[:, :2].copy()
                partial_matrix[:, 2*npeaks+i+1] = partial[:, 2].copy()
            partial_matrix[:, 2*npeaks] = best_heights[-1] * region.indexes if self._baseline else 1
            partial_matrix[:, -1] = 1 + best_params[-1] * region.indexes

        try:
            partial_inverse = np.linalg.solve((partial_matrix.T @ partial_matrix), partial_matrix.T)
        except np.linalg.LinAlgError:
            print('LinAlgError: partial_matrix.T @ partial_matrix is not invertible.')
            print("Use pseudo-inverse instead")
            vals, vecs = np.linalg.eig(partial_matrix.T @ partial_matrix)
            vals, vecs = vals[vals.argsort()[::-1]], vecs[:, vals.argsort()[::-1]]
            inverse = np.zeros((partial_matrix.shape[1], partial_matrix.shape[1]))
            for i in range(partial_matrix.shape[1]):
                if abs(vals[i]) > 1E-4:
                    inverse += (vecs[:, i:i+1] @ vecs[:, i:i+1].T / vals[i])
            partial_inverse = inverse @ partial_matrix.T
        covariance_matrix = np.einsum('ig,gk,jk->ij', partial_inverse, spectrum_covariance_matrix, partial_inverse)
        return covariance_matrix


class EntiretyPeakFitter(Operator):

    def __init__(self, trial: int = 10, label: str | None = None):
        self._trial = trial
        if label is None:
            label = "EntiretyPeakFitter"
        super().__init__(1, label)

    def __run__(self, spectra: list[Spectrum], *args, **kwargs) -> Spectrum:
        fitted = spectra[0].copy()
        fcounts = fitted.copy()

        best_error = 1E8
        guesses = []
        for region in fitted.regions:
            for peak in region.peaks:
                guesses.extend([peak['location'], region.length/region.npeaks/4])
            region.peaks = []

        guesses.append(0)
        npeaks = sum([region.npeaks for region in fitted.regions])
        heights = np.ones((npeaks+1, ))  # mutable object to transfer data without return
        shapes = np.zeros((fitted.shape[0], npeaks+1))

        best_params, best_heights = guesses.copy(), heights.copy()
        for i in range(self._trial):
            params, error = fmin(self._fitgaussian, x0=guesses, args=(fitted.indexes, fitted, shapes, heights, fcounts, npeaks), disp=False, full_output=True)[:2]
            if error <= best_error and heights.min() >= 0:
                best_error = error
                guesses = params * np.random.normal(1, 0.05, len(params))
                best_params = params.copy()
                best_heights = heights.copy()

        for i in range(npeaks):
            center = best_params[i*2]
            stderr = best_params[i*2+1]
            height = best_heights[i]
            for region in fitted.regions:
                if round(center) in region.indexes:
                    region.peaks.append({'location': center, 'height': height,
                                         'stderr': stderr, 'height': height, 'area': height * stderr * (2*np.pi)**0.5})
        for region in fitted.regions:
            region.slope = best_params[-1] * best_heights[-1]
            region.offset = best_heights[-1]
        return fitted

    @staticmethod
    def _fitgaussian(params: np.ndarray, indexes: np.ndarray, counts: np.ndarray, shapes: np.ndarray, heights: np.ndarray, fcounts: np.ndarray, npeaks: int):
        for i in range(npeaks):
            shapes[:, i] = __class__._gaussian(indexes, params[2*i], params[2*i+1])
        shapes[:, -1] = __class__._linear(indexes, params[2*npeaks])
        heights[:] = np.abs(np.abs(np.linalg.lstsq(shapes, counts)[0]))
        fcounts[:] = shapes @ heights
        return np.linalg.norm(fcounts - counts)

    @staticmethod
    def _gaussian(x: np.ndarray, center: float, std: float):
        return np.exp(-(x - center) ** 2 / (2 * std ** 2))

    @staticmethod
    def _linear(x: np.ndarray, slope: float):
        return slope * x + 1


class Deconvolutioner(Operator):

    def __init__(self, trials: int = 10, otrials: int = 10, label: str | None = None):
        self._trials = trials
        self._otrials = otrials
        if label is None:
            label = "Deconvolutioner"
        super().__init__(1, label)

    def __run__(self, spectra: list[Spectrum], *args, **kargs) -> Spectrum:
        deconvoluted = spectra[0].copy()
        for region in deconvoluted.regions:
            deconvoluted[region.indexes] = self._deconvolute(region, spectra[0])
        return deconvoluted

    def _broaden(self, erg_i, erg_j, calib):
        stderr = calib(erg_j) / 2.355
        return np.exp(-(erg_i-erg_j)**2 / (2*stderr**2)) / (2 * np.pi) ** 0.5 / stderr

    def _deconvolute(self, region: Region, spectrum: Spectrum):
        deconvoluted = np.ones(region.length)
        deconvoluted_c = np.ones(region.length) * 2
        deconvoluted_cc = np.ones(region.length) * 3
        ii, jj = region.indexes, region.indexes
        jjmesh, iimesh = np.meshgrid(ii, jj)
        response = self._broaden(iimesh, jjmesh, spectrum.FWHMcal)
        index_o = 0
        while np.abs((deconvoluted_cc+1)/(deconvoluted_c+1)-1).max() > 1E-3 and index_o < self._otrials:
            index_i = 0
            index_o += 1
            deconvoluted_cc = deconvoluted_c.copy()
            deconvoluted_c = deconvoluted_c ** 1.4
            while np.abs((deconvoluted+1)/(deconvoluted_c+1)-1).max() > 1E-3 and index_i < self._trials:
                index_i += 1
                deconvoluted = deconvoluted_c.copy()
                deconvoluted_c = (response.T @ spectrum[region.indexes]) * deconvoluted / ((response.T @ response) @ deconvoluted)

        return deconvoluted_c
