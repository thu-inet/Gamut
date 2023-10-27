import numpy as np
from scipy.optimize import fmin, minimize, Bounds
from copy import deepcopy
from time import strftime, ctime

from .Spectrum import Spectrum
from .Operator import Operator
from .PeakRegion import Region


class AverageClassicPeakAreaCalculator(Operator):

    def __init__(self, times: int, label: str = None):
        """
        The most basic PeakAreaCalculator.

        :param base_mode: mode to strip counts of baseline.
        :param label: label for the operator
        """
        self._times = times
        if label is None:
            label = f"FPAreaCalculator[T{self._times}]"
        super().__init__(1, label)

    def __run__(self, spectra: list[Spectrum], *args, **kargs) -> Spectrum:
        calculated = spectra[0].copy()
        for region in calculated.regions:
            self._calculate(region, calculated)
        return calculated

    def _split_regions(self, region: Region, spectrum: Spectrum):
        splited_region_indexes = [(spectrum[region.peaks[i]['location']: region.peaks[i+1]['location']]).argmin()+region.peaks[i]['location'] for i in range(region.npeaks-1)]
        splited_region_indexes = [region.left] + splited_region_indexes + [region.right+1]
        return splited_region_indexes

    def _calculate(self, region: Region, spectrum: Spectrum):
        if region.npeaks == 0:
            raise ValueError("Region has no peaks")
        elif region.npeaks == 1:
            peak_areas = []
            for i in range(self._times):
                left = max(region.left - i, 0)
                right = min(region.right + i, len(spectrum)-1)
                total = spectrum[left: right+1].sum()
                baseline = (spectrum[left] + spectrum[right]) * (right-left+1) / 2
                peak_areas.append(total - baseline)
            region.peaks[0]['FParea'] = sum(peak_areas) / len(peak_areas)
        else:
            splited_region_indexes = self._split_regions(region, spectrum)
            for i, peak in enumerate(region.peaks):
                total = spectrum[splited_region_indexes[i]: splited_region_indexes[i+1]+1].sum()
                baseline = (spectrum[splited_region_indexes[i]] + spectrum[min(splited_region_indexes[i+1], len(spectrum)-1)]) * (splited_region_indexes[i+1] - splited_region_indexes[i] + 1) / 2
                # print(total, baseline, spectrum[splited_region_indexes[i]], spectrum[splited_region_indexes[i+1]])
                peak['FParea'] = total - baseline


# class AverageCovellPeakAreaCalculator(Operator):

#     def __init__(self, init_width: int, times: int, label: str = None):
#         self._init_width = init_width
#         self._times = times
#         if label is None:
#             label = f"CovellCalculator[T{self._times}]"
#         super().__init__(1, label)

#     def __run__(self, spectra: list[Spectrum], *args, **kargs) -> Spectrum:
#         calculated = spectra[0].copy()
#         for region in calculated.regions:
#             peak_areas = self._calculate(region)
#             region.area = sum(peak_areas) / len(peak_areas)
#         return calculated

#     def _calculate(self, region: Region):
#         spectrum = region.spectrum
#         peak_areas = []
#         for width in range(self._init_width, self._init_width+self._times+1):
#             total = spectrum[region.median('highest') - width: region.median('highest') + width + 1].sum()
#             base = (spectrum[region.median('highest') - width] + spectrum[region.median('highest') + width]) * (2 * width + 1) / 2
#             peak_areas.append(total - base)
#         return peak_areas


# class AverageWassonPeakAreaCalculator(Operator):

#     def __init__(self, init_width: int, times: int, label: str = None):
#         self._init_width = init_width
#         self._times = times
#         if label is None:
#             label = f"CovellCalculator[T{self._times}]"
#         super().__init__(1, label)

#     def __run__(self, spectra: list[Spectrum], *args, **kargs) -> Spectrum:
#         calculated = spectra[0].copy()
#         for peak in calculated.regions:
#             peak_areas = self._calculate(peak)
#             if len(peak_areas) > 0:
#                 peak.area = sum(peak_areas) / len(peak_areas)
#             else:
#                 peak.area = 0
#         return calculated

#     def _calculate(self, peak: Region):
#         spectrum = peak.spectrum
#         peak_areas = []
#         baseline = np.zeros(spectrum.shape)
#         baseline[peak.indexes] = np.linspace(spectrum[peak.left], spectrum[peak.right], peak.length)
#         for width in range(self._init_width, self._init_width+self._times+1):
#             if (peak.median('highest') - width < peak.left) or (peak.median('highest') + width > peak.right):
#                 break
#             total = spectrum[peak.median('highest') - width: peak.median('highest') + width + 1].sum()
#             base = (baseline[peak.median('highest') - width] + baseline[peak.median('highest') + width]) * (2 * width + 1) / 2
#             peak_areas.append(total - base)
#         return peak_areas


class RegionPeakFitter(Operator):

    def __init__(self, trial: int = 10, equal_width: bool = True,
                 baseline: bool = True, variable_npeaks: bool = False, label: str = None):
        self._trial = trial
        self._baseline = baseline
        self._equal_width = equal_width
        self._variable_npeaks = variable_npeaks
        if label is None:
            label = "RegionPeakFitter"
        super().__init__(1, label)

    def __run__(self, spectra: list[Spectrum], *args, **kwargs) -> Spectrum:
        fitted = spectra[0].copy()
        for region in fitted.regions:
            print(f"Fitting Region: {region.left}~{region.right}, NPeaks={region.npeaks}, time={strftime(ctime())}")
            best_params, best_heights, _, _, npeaks = self._fit(region, fitted)
            if self._equal_width:
                region.peaks = []
                for i in range(len(best_heights)-1):
                    if round(best_params[i]) in region.indexes:
                        peak = {}
                        peak['center'] = best_params[i]
                        peak['location'] = round(best_params[i])
                        peak['height'] = best_heights[i]
                        peak['stderr'] = best_params[-2]
                        peak['area'] = peak['height'] * peak['stderr'] * (2*np.pi)**0.5
                        region.peaks.append(peak)
                region.slope = best_params[-1] * best_heights[-1]
                region.offset = best_heights[-1]
            else:
                region.peaks = []
                for i in range(len(best_heights)-1):
                    if round(best_params[2*i]) in region.indexes:
                        peak = {}
                        peak['center'] = best_params[2*i]
                        peak['location'] = round(best_params[2*i])
                        peak['height'] = best_heights[i]
                        peak['stderr'] = best_params[2*i+1]
                        peak['area'] = peak['height'] * peak['stderr'] * (2*np.pi)**0.5
                        region.peaks.append(peak)
                region.slope = best_params[-1] * best_heights[-1]
                region.offset = best_heights[-1]
                region.peaks = sorted(region.peaks, key=lambda k: k['location'])
        return fitted

    def _fit(self, region: Region, fitted: Spectrum):
        best_error = 1E8
        npeaks = len(region.peaks) + 1
        heights = np.ones(npeaks)  # mutable object to transfer data without return
        shapes = np.zeros((region.length, npeaks))
        fcounts = fitted[region.indexes].copy()

        if self._equal_width:
            guesses = np.array([peak['location'] for peak in region.peaks])
            guesses = np.hstack((guesses, [region.length/npeaks/4, 0]))
        else:
            guesses = np.array([[peak['location'], region.length/npeaks/4] for peak in region.peaks]).ravel()
            guesses = np.hstack((guesses, 0))

        for i in range(self._trial):
            params, error = fmin(self._fitgaussian, x0=guesses,
                                 args=(region.indexes, fitted[region.indexes], shapes, heights, fcounts, npeaks),
                                 disp=False, full_output=True)[:2]
            if (error <= best_error and heights[:-1].min() >= 0):
                best_error = error
                guesses = params * np.random.normal(1, 0.05, len(params))
                best_params = params.copy()
                best_heights = heights.copy()
                best_shapes = shapes.copy()
                best_fcounts = fcounts.copy()
        return best_params, best_heights, best_shapes, best_fcounts, npeaks

    # def error_estimate(self, region: 'Region', fitted: np.ndarray, npeaks: int):
    #     """
    #     Estimate the error of the fitting.

    #     Parameters
    #     ----------
    #     region : Region
    #         The region to be fitted.
    #     fitted : np.ndarray
    #         The fitted data.
    #     npeaks : int
    #         The number of peaks.

    #     Returns
    #     -------
    #     float
    #         The error of the fitting.
    #     """
    #     return np.sum((fitted - region.counts) ** 2) / (region.counts.size - npeaks)

    def _fitgaussian(self, params: np.ndarray, indexes: np.ndarray, counts: np.ndarray,
                     shapes: np.ndarray, heights: np.ndarray, fcounts: np.ndarray, npeaks: int):
        if self._equal_width:
            for i in range(npeaks-1):
                shapes[:, i] = __class__._gaussian(indexes, params[i], params[-2])
        else:
            for i in range(npeaks-1):
                shapes[:, i] = __class__._gaussian(indexes, params[2*i], params[2*i+1])
        shapes[:, -1] = __class__._linear(indexes, params[-1]) * self._baseline
        heights[:] = np.linalg.lstsq(shapes, counts)[0]
        heights[:-1] = np.abs(heights[:-1])
        fcounts[:] = shapes @ heights
        return np.linalg.norm(fcounts - counts, ord=2)

    @staticmethod
    def _gaussian(x: np.ndarray, center: float, std: float):
        return np.exp(-(x - center) ** 2 / (2 * std ** 2))

    @staticmethod
    def _linear(x: np.ndarray, slope: float):
        return slope * x + 1


# class EntiretyPeakFitter(Operator):

#     def __init__(self, trial: int = 10, label: str = None):
#         self._trial = trial
#         if label is None:
#             label = "EntiretyPeakFitter"
#         super().__init__(1, label)

#     def __run__(self, spectra: list[Spectrum], *args, **kwargs) -> Spectrum:
#         fitted = spectra[0].copy()
#         best_error = 1000
#         npeaks = len(fitted.regions)
#         heights = np.ones((npeaks+1, ))  # mutable object to transfer data without return
#         shapes = np.zeros((fitted.shape[0], npeaks+1))
#         fcounts = fitted.copy()
#         guess_params = []
#         for peak in fitted.regions:
#             guess_params.append(peak.median())
#             guess_params.append(peak.median()-peak.left)
#         guess_params.append(0)
#         guess_params = np.array(guess_params, dtype=np.float64).ravel()
#         for _ in range(self._trial):
#             params, error = fmin(self._fitgaussian, x0=guess_params, args=(fitted.indexes, fitted, shapes, heights, fcounts, npeaks), disp=False, full_output=True)[:2]
#             if error <= best_error and heights.min() >= 0:
#                 best_error = error
#                 guess_params = params * np.random.normal(1, 0.05, len(params))
#                 best_params = params.copy()
#                 best_shapes = shapes.copy()
#                 best_heights = heights.copy()
#                 best_fcounts = fcounts.copy()
#         fitted[:] = fcounts
#         fitted.heights = best_heights
#         fitted.peak_params = best_params
#         fitted.npeaks = npeaks
#         return fitted

#     @staticmethod
#     def _fitgaussian(params: np.ndarray, indexes: np.ndarray, counts: np.ndarray, shapes: np.ndarray, heights: np.ndarray, fcounts: np.ndarray, npeaks: int):
#         for i in range(npeaks):
#             shapes[:, i] = __class__._gaussian(indexes, params[2*i], params[2*i+1])
#         shapes[:, -1] = __class__._linear(indexes, params[2*npeaks])
#         heights[:] = np.abs(np.abs(np.linalg.lstsq(shapes, counts)[0]))
#         fcounts[:] = shapes @ heights
#         return np.linalg.norm(fcounts - counts)

#     @staticmethod
#     def _gaussian(x: np.ndarray, center: float, std: float):
#         return np.exp(-(x - center) ** 2 / (2 * std ** 2))

#     @staticmethod
#     def _linear(x: np.ndarray, slope: float):
#         return slope * x + 1


class Deconvolutioner(Operator):

    def __init__(self, trials: int = 10, otrials: int = 10, label: str = None):
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
            # print(index_o, deconvoluted_c.sum())
            # deconvoluted_c = deconvoluted_c / deconvoluted_c.sum() * spectrum[region.indexes].sum()
            # print(sum(deconvoluted), sum(deconvoluted_c), sum(deconvoluted_cc))
        return deconvoluted_c


if __name__ == "__main__":

    import matplotlib.pyplot as plt
    # from Spectrum import simuspecs
    # from Smoother import SavitzkySmoother
    # from PeakSearcher import GaussPeakSearcher, DifferentialPeakSearcher
    # from OtherOperator import IterativeNonEquilibriumWhittakerStripper
    # from Operator import Pipe, PipeNet

    # sav = SavitzkySmoother(2, 5)
    # strp = IterativeNonEquilibriumWhittakerStripper(1E4, 40)
    # gauss = GaussPeakSearcher(2, 0.05)
    # diff = DifferentialPeakSearcher(2, 3, 2, 4, 0.5)
    # # pipe = Pipe([sav, gauss])

    # simu0 = simuspecs['synthesized'].copy()
    # # simu0 = simuspecs['doublepeak_normal_narrow'].copy()
    # stripped = gauss(sav(simu0) - strp(sav(simu0)))
    # simu = simu0.copy()
    # simu.regions = stripped.regions
    # # simu[:] = sav(simu0) - strp(simu)
    # simu[:] =  simu - strp(simu0)

    # # simu.plot()
    # # for spec in pipe.list_spectrum:
    # #     spec.plot()
    # # plt.legend()
    # # plt.show()

    # # fpa = AverageClassicPeakAreaCalculator(5)
    # # covell = AverageCovellPeakAreaCalculator(24, 10)
    # # wasson = AverageWassonPeakAreaCalculator(15, 10)
    # # area1 = fpa(simu)
    # # area2 = covell(simu)
    # # area3 = wasson(simu)

    # # for area in [area1, area2, area3]:
    # #     print(f"Area by {area.label} is", area.regions[0].area)

    # fiter = SinglePeakFitter()
    # fit = fiter(simu)

    # # for peak in fit.regions:
    # #     print(peak.fit_results[0][0], peak.fit_results[0][1], peak.fit_results[0][1] * peak.fit_results[2][0] * (np.pi * 2)**0.5)

    # simu0.plot()
    # base = strp(simu0)
    # base.plot()
    # bbase = base.copy()
    # for peak in fit.regions:
    #     if round(peak.fit_results[0][0]) in peak.indexes and peak.fit_results[2][0] <= 2*simu[peak.indexes].max():
    #         indexes = np.arange(int(peak.fit_results[0][0]-5*peak.fit_results[0][1]), min(int(peak.fit_results[0][0]+5*peak.fit_results[0][1]), simu.length-1))
    #         shapes = SinglePeakFitter._gaussian(indexes, peak.fit_results[0][0], peak.fit_results[0][1]) * peak.fit_results[2][0]
    #         plt.fill_between(indexes, base[indexes], base[indexes]+shapes, alpha=0.5)
    #         bbase[indexes] += shapes
    # bbase.plot(linestyle='--')
    # plt.show()


    # # fit.plot()
    # # fit.plot_peaks()
    # # for peak in fit.regions:
    # #     if round(peak.fit_results[0][0]) in peak.indexes and peak.fit_results[2][0] <= 2*simu[peak.indexes].max():
    # #         indexes = np.arange(int(peak.fit_results[0][0]-5*peak.fit_results[0][1]), int(peak.fit_results[0][0]+5*peak.fit_results[0][1]))
    # #         shapes = SinglePeakFitter._gaussian(indexes, peak.fit_results[0][0], peak.fit_results[0][1]) * peak.fit_results[2][0]
    # #         plt.fill_between(indexes, np.zeros(indexes.shape), shapes, alpha=0.5)
    # #         # plt.plot(peak.indexes, peak.fit_results[-1])
    # #         # plt.fill_between(peak.indexes, np.zeros(peak.length), peak.fit_results[1][:, 0] * peak.fit_results[2][0])
    # #         # plt.plot(peak.indexes, peak.fit_results[1][:, 1] * peak.fit_results[2][1])
    # # plt.show()


    # # enfit = EntiretyPeakFitter()
    # # fit = enfit(simu) 
    # # fit += strp(simu0)
    # # simu0.plot()
    # # fit.plot()
    # # for i in range(fit.npeaks):
    # #     plt.plot(enfit._gaussian(fit.indexes, fit.peak_params[2*i], fit.peak_params[2*i+1])*fit.heights[i] + strp(simu0))
    # # plt.show()