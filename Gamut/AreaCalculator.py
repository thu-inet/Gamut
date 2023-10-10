import numpy as np
from scipy.optimize import fmin
from copy import deepcopy

from Spectrum import Spectrum
from Operator import Operator
from PeakRegion import PeakRegion


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
        for peak in calculated.peaks:
            peak_areas = self._calculate(peak)
            peak.area = sum(peak_areas) / len(peak_areas)
        return calculated

    def _calculate(self, peak: PeakRegion):
        spectrum = peak.spectrum
        peak_areas = []
        for i in range(self._times):
            total = spectrum[peak.left-i: peak.right+1+i].sum()
            base = (spectrum[peak.left-i] + spectrum[peak.right+i]) * (peak.length + 2*i) / 2
            peak_areas.append(total - base)
            # print((peak.length + 2*i), total, base, peak_areas[-1])
        return peak_areas


class AverageCovellPeakAreaCalculator(Operator):

    def __init__(self, init_width: int, times: int, label: str = None):
        self._init_width = init_width
        self._times = times
        if label is None:
            label = f"CovellCalculator[T{self._times}]"
        super().__init__(1, label)

    def __run__(self, spectra: list[Spectrum], *args, **kargs) -> Spectrum:
        calculated = spectra[0].copy()
        for peak in calculated.peaks:
            peak_areas = self._calculate(peak)
            peak.area = sum(peak_areas) / len(peak_areas)
        return calculated

    def _calculate(self, peak: PeakRegion):
        spectrum = peak.spectrum
        peak_areas = []
        for width in range(self._init_width, self._init_width+self._times+1):
            total = spectrum[peak.median('highest') - width: peak.median('highest') + width + 1].sum()
            base = (spectrum[peak.median('highest') - width] + spectrum[peak.median('highest') + width]) * (2 * width + 1) / 2
            peak_areas.append(total - base)
        return peak_areas


class AverageWassonPeakAreaCalculator(Operator):

    def __init__(self, init_width: int, times: int, label: str = None):
        self._init_width = init_width
        self._times = times
        if label is None:
            label = f"CovellCalculator[T{self._times}]"
        super().__init__(1, label)

    def __run__(self, spectra: list[Spectrum], *args, **kargs) -> Spectrum:
        calculated = spectra[0].copy()
        for peak in calculated.peaks:
            peak_areas = self._calculate(peak)
            if len(peak_areas) > 0:
                peak.area = sum(peak_areas) / len(peak_areas)
            else:
                peak.area = 0
        return calculated

    def _calculate(self, peak: PeakRegion):
        spectrum = peak.spectrum
        peak_areas = []
        baseline = np.zeros(spectrum.shape)
        baseline[peak.indexes] = np.linspace(spectrum[peak.left], spectrum[peak.right], peak.length)
        for width in range(self._init_width, self._init_width+self._times+1):
            if (peak.median('highest') - width < peak.left) or (peak.median('highest') + width > peak.right):
                break
            total = spectrum[peak.median('highest') - width: peak.median('highest') + width + 1].sum()
            base = (baseline[peak.median('highest') - width] + baseline[peak.median('highest') + width]) * (2 * width + 1) / 2
            peak_areas.append(total - base)
        return peak_areas


class SinglePeakFitter(Operator):

    def __init__(self, trial: int = 10, baseline: str = "None", label: str = None):
        self._trial = trial
        self._baseline = baseline
        if label is None:
            label = "SinglePeakFitter"
        super().__init__(1, label)

    def __run__(self, spectra: list[Spectrum], *args, **kwargs) -> Spectrum:
        fitted = spectra[0].copy()
        for peak in fitted.peaks:
            fit_results = self._fit(peak, fitted)
            peak.fit_results = deepcopy(fit_results)
        return fitted

    def _fit(self, peak: PeakRegion, fitted: Spectrum):
        best_error = 1000
        heights = np.ones((2, ))  # mutable object to transfer data without return
        shapes = np.zeros((fitted[peak.indexes].shape[0], 2))
        fcounts = fitted[peak.indexes].copy()
        guess_params = np.array([peak.median(), peak.median()-peak.left, 0], dtype=np.float64).ravel()
        for _ in range(self._trial):
            params, error = fmin(self._fitgaussian, x0=guess_params, args=(peak.indexes, fitted[peak.indexes], shapes, heights, fcounts), disp=False, full_output=True)[:2]
            if error <= best_error and heights.min() >= 0:
                best_error = error
                guess_params = params * np.random.normal(1, 0.05, len(params))
                best_params = params.copy()
                best_shapes = shapes.copy()
                best_heights = heights.copy()
                best_fcounts = fcounts.copy()
        return best_params, best_shapes, best_heights, best_fcounts

    @staticmethod
    def _fitgaussian(params: np.ndarray, indexes: np.ndarray, counts: np.ndarray, shapes: np.ndarray, heights: np.ndarray, fcounts: np.ndarray):
        shapes[:, 0] = __class__._gaussian(indexes, params[0], params[1])
        shapes[:, 1] = __class__._linear(indexes, params[2]) * 0
        heights[:] = np.abs(np.abs(np.linalg.lstsq(shapes, counts)[0]))
        fcounts[:] = shapes @ heights
        return np.linalg.norm(fcounts - counts)

    @staticmethod
    def _gaussian(x: np.ndarray, center: float, std: float):
        return np.exp(-(x - center) ** 2 / (2 * std ** 2))

    @staticmethod
    def _linear(x: np.ndarray, slope: float):
        return slope * x + 1


class EntiretyPeakFitter(Operator):

    def __init__(self, trial: int = 10, label: str = None):
        self._trial = trial
        if label is None:
            label = "EntiretyPeakFitter"
        super().__init__(1, label)

    def __run__(self, spectra: list[Spectrum], *args, **kwargs) -> Spectrum:
        fitted = spectra[0].copy()
        best_error = 1000
        npeaks = len(fitted.peaks)
        heights = np.ones((npeaks+1, ))  # mutable object to transfer data without return
        shapes = np.zeros((fitted.shape[0], npeaks+1))
        fcounts = fitted.copy()
        guess_params = []
        for peak in fitted.peaks:
            guess_params.append(peak.median())
            guess_params.append(peak.median()-peak.left)
        guess_params.append(0)
        guess_params = np.array(guess_params, dtype=np.float64).ravel()
        for _ in range(self._trial):
            params, error = fmin(self._fitgaussian, x0=guess_params, args=(fitted.indexes, fitted, shapes, heights, fcounts, npeaks), disp=False, full_output=True)[:2]
            if error <= best_error and heights.min() >= 0:
                best_error = error
                guess_params = params * np.random.normal(1, 0.05, len(params))
                best_params = params.copy()
                best_shapes = shapes.copy()
                best_heights = heights.copy()
                best_fcounts = fcounts.copy()
        # print(best_params)
        fitted[:] = fcounts
        fitted.heights = best_heights
        fitted.peak_params = best_params
        fitted.npeaks = npeaks
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
    

if __name__ == "__main__":

    import matplotlib.pyplot as plt
    from Spectrum import simuspecs
    from Smoother import SavitzkySmoother
    from PeakSearcher import GaussPeakSearcher, DifferentialPeakSearcher
    from OtherOperator import IterativeNonEquilibriumWhittakerStripper
    from Operator import Pipe, PipeNet

    sav = SavitzkySmoother(2, 5)
    strp = IterativeNonEquilibriumWhittakerStripper(1E4, 40)
    gauss = GaussPeakSearcher(2, 0.05)
    diff = DifferentialPeakSearcher(2, 3, 2, 4, 0.5)
    # pipe = Pipe([sav, gauss])

    simu0 = simuspecs['synthesized'].copy()
    # simu0 = simuspecs['doublepeak_normal_narrow'].copy()
    stripped = gauss(sav(simu0) - strp(sav(simu0)))
    simu = simu0.copy()
    simu.peaks = stripped.peaks
    # simu[:] = sav(simu0) - strp(simu)
    simu[:] =  simu - strp(simu0)

    # simu.plot()
    # for spec in pipe.list_spectrum:
    #     spec.plot()
    # plt.legend()
    # plt.show()

    # fpa = AverageClassicPeakAreaCalculator(5)
    # covell = AverageCovellPeakAreaCalculator(24, 10)
    # wasson = AverageWassonPeakAreaCalculator(15, 10)
    # area1 = fpa(simu)
    # area2 = covell(simu)
    # area3 = wasson(simu)

    # for area in [area1, area2, area3]:
    #     print(f"Area by {area.label} is", area.peaks[0].area)

    fiter = SinglePeakFitter()
    fit = fiter(simu)
    
    # for peak in fit.peaks:
    #     print(peak.fit_results[0][0], peak.fit_results[0][1], peak.fit_results[0][1] * peak.fit_results[2][0] * (np.pi * 2)**0.5)

    simu0.plot()
    base = strp(simu0)
    base.plot()
    bbase = base.copy()
    for peak in fit.peaks:
        if round(peak.fit_results[0][0]) in peak.indexes and peak.fit_results[2][0] <= 2*simu[peak.indexes].max():
            indexes = np.arange(int(peak.fit_results[0][0]-5*peak.fit_results[0][1]), min(int(peak.fit_results[0][0]+5*peak.fit_results[0][1]), simu.length-1))
            shapes = SinglePeakFitter._gaussian(indexes, peak.fit_results[0][0], peak.fit_results[0][1]) * peak.fit_results[2][0]
            plt.fill_between(indexes, base[indexes], base[indexes]+shapes, alpha=0.5)
            bbase[indexes] += shapes
    bbase.plot(linestyle='--')
    plt.show()


    # fit.plot()
    # fit.plot_peaks()
    # for peak in fit.peaks:
    #     if round(peak.fit_results[0][0]) in peak.indexes and peak.fit_results[2][0] <= 2*simu[peak.indexes].max():
    #         indexes = np.arange(int(peak.fit_results[0][0]-5*peak.fit_results[0][1]), int(peak.fit_results[0][0]+5*peak.fit_results[0][1]))
    #         shapes = SinglePeakFitter._gaussian(indexes, peak.fit_results[0][0], peak.fit_results[0][1]) * peak.fit_results[2][0]
    #         plt.fill_between(indexes, np.zeros(indexes.shape), shapes, alpha=0.5)
    #         # plt.plot(peak.indexes, peak.fit_results[-1])
    #         # plt.fill_between(peak.indexes, np.zeros(peak.length), peak.fit_results[1][:, 0] * peak.fit_results[2][0])
    #         # plt.plot(peak.indexes, peak.fit_results[1][:, 1] * peak.fit_results[2][1])
    # plt.show()


    # enfit = EntiretyPeakFitter()
    # fit = enfit(simu) 
    # fit += strp(simu0)
    # simu0.plot()
    # fit.plot()
    # for i in range(fit.npeaks):
    #     plt.plot(enfit._gaussian(fit.indexes, fit.peak_params[2*i], fit.peak_params[2*i+1])*fit.heights[i] + strp(simu0))
    # plt.show()