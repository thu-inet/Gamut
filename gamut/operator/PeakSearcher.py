import numpy as np

from typing import Literal, Callable
from multiprocessing import Pool
from scipy.signal import argrelextrema

from ..utils import Differential
from ..classes import Region

from .Operator import Operator
from ..spectrum.Spectrum import Spectrum
from .Smoother import CentroidSmoother

class PeakSearcher(Operator):

    def __init__(self,
                 shift: int = 0,
                 min_length: int = 4,
                 min_height: float = 30,
                 min_height_ratio: float = 2,
                 min_area: float = 30,
                 min_area_ratio: float = 1,
                 merge_peaks: bool = True,
                 merge_mode: Literal["No", "Old", "New"] | None = None,
                 print_deprecate_info: bool = False,
                 note_deprecate_info: bool = False,
                 label: str | None = None):
        self._shift = shift
        self._min_length = min_length
        self._min_height = min_height
        self._min_height_ratio = min_height_ratio
        self._min_area = min_area
        self._min_area_ratio = min_area_ratio
        if merge_mode is None:
            merge_mode = "New"
        self._merge_peaks = merge_peaks
        self._merge_mode = merge_mode
        self._print_deprecate_info = print_deprecate_info
        self._note_deprecate_info = note_deprecate_info
        super().__init__(1, label)

    @staticmethod
    def _split(data: np.ndarray) -> list[Region]:
        """
        Split a bool sequence into pieces of consecutive True.
        """
        regions = []

        start, i = None, 0
        for i, d in enumerate(data):
            if d:
                if start is None:
                    start = i
            else:
                if start is not None:
                    regions.append(Region(start, i-1))
                    start = None
        if start is not None:
            regions.append(Region(start, i))
        return regions

    @staticmethod
    def _expand(regions: list[Region], spectrum: Spectrum, mode: Literal["max_area", "min_height"]="min_height") -> list[Region]:
        smooothed = CentroidSmoother(4)(spectrum)

        if mode == "max_area":
            for region in regions:
                while True:
                    total, baseline = spectrum.estimate_area(region)
                    area = total - baseline
                    expanded_region = [Region(max(region.left-1, 0), min(region.right+0, spectrum.length-1)),
                                    Region(max(region.left-0, 0), min(region.right+1, spectrum.length-1)),
                                    Region(max(region.left-1, 0), min(region.right+1, spectrum.length-1)),
                                    Region(max(region.left-2, 0), min(region.right+1, spectrum.length-1)),
                                    Region(max(region.left-1, 0), min(region.right+2, spectrum.length-1)),
                                    Region(max(region.left-2, 0), min(region.right+2, spectrum.length-1))]
                    expanded_area = [spectrum.estimate_area(region2)[0]-spectrum.estimate_area(region2)[1] for region2 in expanded_region]
                    expanded_area = np.array(expanded_area)
                    if np.any(expanded_area > area):
                        region2 = expanded_region[np.argmax(expanded_area)]
                        region.left = region2.left
                        region.right = region2.right
                    else:
                        break
        elif mode == "min_height":
            for region in regions:
                if len(argrelextrema(smooothed[ :region.left+2], np.less)[0]) > 0:
                    region.left = argrelextrema(smooothed[ :region.left+2], np.less)[0][-1]
                if len(argrelextrema(smooothed[region.right-2: ], np.less)[0]) > 0:
                    region.right = argrelextrema(smooothed[region.right-2: ], np.less)[0][0] + region.right
        else:
            raise ValueError(f"mode {mode} not supported.")

        return regions

    def _merge(self, regions: list[Region], searched: Spectrum, merge_mode: Literal["No", "Old", "New"] | None = None) -> list[Region]:

        if not self._merge_peaks:
            return regions

        if merge_mode is not None:
            merge_mode = merge_mode
        else:
            merge_mode = self._merge_mode

        if merge_mode == "No":
            regions = regions
        elif merge_mode == "Old":
            regions = searched.regions + regions
        elif merge_mode == "New":
            regions = regions + searched.regions
        else:
            raise ValueError(f"mode {merge_mode} not supported.")
        
        if len(regions) == 0:
            return regions

        max_index = max([region.right for region in regions])
        united_indexes = np.zeros(max_index+1)
        united_peaks = []
        for region in regions:
            united_indexes[region.indexes] = 1
            united_peaks.extend(region.peaks)
        merged_regions = self._split(united_indexes)
        # united_peaks = sorted(united_peaks, key=lambda x: x['location'])
        
        # remove repeated peaks
        new_peaks = []
        peak_locations = []
        for peak in united_peaks:
            if peak['location'] not in peak_locations:
                new_peaks.append(peak)
                peak_locations.extend([i for i in range(peak['location']-1, peak['location']+2)])

        for region in merged_regions:
            region.peaks = [peak for peak in new_peaks if (peak['location'] <= region.right) and (peak['location'] >= region.left)]
        return merged_regions

    def _correct(self, regions: list[Region], spectrum: Spectrum, shift: int | None = None) -> list[Region]:
        if shift is None:
            shift = self._shift
        regions_copy = []
        for peak in regions:
            peak_copy = peak.copy()
            peak_copy.left = max(peak.left+shift, 0)
            peak_copy.right = min(peak.right+shift, spectrum.shape[0] - 1)
            regions_copy.append(peak_copy)
        return regions_copy

    def generic_filter(self, func: Callable, info: str):
        def remove(self, spectrum: Spectrum, regions: list[Region]) -> list[Region]:
            regions_copy = []
            for region in regions:
                if func(spectrum, region):
                    regions_copy.append(region)
                    if self._print_deprecate_info:
                        print(f"PeakSearcher: region at {region.left:>5d}~{region.right:<5d} (Energy={spectrum.ergcal(region.left):>8.2f}~{spectrum.ergcal(region.right):<8.2f} keV) passed in test for {info}")
                else:
                    if self._print_deprecate_info:
                        print(f"PeakSearcher: region at {region.left:>5d}~{region.right:<5d} (Energy={spectrum.ergcal(region.left):>8.2f}~{spectrum.ergcal(region.right):<8.2f} keV) failed in test for {info}")
            return regions_copy
        return remove

    def _filter_by_length(self, regions: list[Region], spectrum: Spectrum,
                       min_length: int | None = None) -> list[Region]:
        if min_length is None:
            min_length = self._min_length
        func = lambda spectrum, region: region.length >= min_length
        regions = self.generic_filter(func, f"too short length")(self, spectrum, regions)
        return regions

    def _filter_by_height(self, regions: list[Region], spectrum: Spectrum,
                       min_height: float | None = None, min_height_ratio: float | None = None) -> list[Region]:
        if min_height is None:
            min_height = self._min_height
        func = lambda spectrum, region: spectrum[region.indexes].max() >= min_height + spectrum[region.indexes].min()
        regions = self.generic_filter(func, f"too low peak height")(self, spectrum, regions)

        if min_height_ratio is None:
            min_height_ratio = self._min_height_ratio
        func = lambda spectrum, region: (spectrum[region.indexes].max() - spectrum[region.indexes].min()) >= min_height_ratio * (spectrum[region.indexes].min() + 1)**0.5
        regions = self.generic_filter(func, f"too low peak height ratio")(self, spectrum, regions)

        return regions

    def _filter_by_area(self, regions: list[Region], spectrum: Spectrum,
                     min_area: float | None = None, min_area_ratio: float | None = None) -> list[Region]:
        if min_area is None:
            min_area = self._min_area
        regions = self.generic_filter(lambda spectrum, region: spectrum.estimate_area(region)[0] - spectrum.estimate_area(region)[1] >= min_area, f"too low peak area")(self, spectrum, regions)
        
        if min_area_ratio is None:
            min_area_ratio = self._min_area_ratio
        regions = self.generic_filter(lambda spectrum, region: (spectrum.estimate_area(region)[0] - spectrum.estimate_area(region)[1]) / (spectrum.estimate_area(region)[0] + 1)**0.5 >= min_area_ratio, f"too low peak area ratio")(self, spectrum, regions)
        return regions
    

class GaussPeakSearcher(PeakSearcher):
    '''
    if even o, o = o/2, g(i) = c(i+o) * c(i-o) / c(i+o+2) / c(i-o-2)
    if odd o, o = (o-1)/2 g(i) = c(i+o+1) * c(i-o) / c(i+o+3) / c(i-o-2)

    '''
    def __init__(self, order: int, threshold: float, label: str | None = None):
        self._order = order
        self._threshold = threshold
        if label is None:
            label = f'GaussPeakSearcher[O{self._order}]'
        super().__init__(shift=round(self._order/2-1), label=label)

    def __run__(self, spectra: list[Spectrum], *args, **kargs) -> Spectrum:
        searched = spectra[0].copy()
        padded = self._transform(searched)
        regions = self._split(padded > self._threshold)
        regions = self._correct(regions, searched)
        for region in regions:
            region.peaks.append({'location': int((region.left+region.right) / 2)})
        regions = self._filter_by_length(regions, searched)
        regions = self._expand(regions, searched)
        regions = self._merge(regions, searched)
        searched.regions = regions
        return searched

    def _transform(self, searched: Spectrum) -> np.ndarray:
        padded = np.maximum(np.pad(searched, (2, self._order), mode='constant'), 1)
        gaussed = padded[2: -self._order] / padded[0: -self._order-2] / padded[self._order+2:] * padded[self._order: -2]
        rectified = (gaussed - 1) * searched ** 0.5
        return rectified


class DifferentialPeakSearcher(PeakSearcher):

    def __init__(self, poly_order: int = 1, hwidth: int = 1, derive_order: int = 1, zero_width: int = 4, label: str | None = None):
        self._poly_order = poly_order
        self._hwidth = hwidth
        self._derive_order = derive_order
        self._zero_width = zero_width
        if label is None:
            label = f'DifferentialPeakSearcher[O{self._derive_order}]'
        super().__init__(label=label)

    def __run__(self, spectra: list[Spectrum], *args, **kargs) -> Spectrum:
        searched = spectra[0].copy()
        differentiated = self._transform(searched)
        zeros = np.zeros_like(differentiated, dtype=bool)
        for i in range(self._zero_width, differentiated.shape[0]-self._zero_width):
            if self._derive_order % 2 == 1:
                zeros[i] = np.all(differentiated[i-self._zero_width-1: i-1] >= differentiated[i]) and np.all(differentiated[i+1: i+self._zero_width+1] <= differentiated[i])
            elif self._derive_order % 2 == 0:
                zeros[i] = np.all(differentiated[i-self._zero_width-1: i-1] >= differentiated[i]) and np.all(differentiated[i+1: i+self._zero_width+1] >= differentiated[i])
        regions = self._split(zeros)
        for region in regions:
            region.peaks.append({'location': round((region.left+region.right) / 2)})
            region.left = max(region.left - self._zero_width - 1, 0)
            region.right = min(region.right + self._zero_width + 1, searched.length - 1)
        regions = self._filter_by_length(regions, searched)
        regions = self._filter_by_height(regions, searched)
        regions = self._expand(regions, searched)
        regions = self._merge(regions, searched)
        regions = self._filter_by_area(regions, searched)
        searched.regions = regions
        return searched

    def _transform(self, searched: Spectrum) -> np.ndarray:
        padded = np.pad(searched, (self._hwidth, self._hwidth), mode='edge')
        differentiated, _ = Differential(padded, self._poly_order, self._hwidth, self._derive_order)
        sliced = differentiated[self._hwidth: -self._hwidth]
        return sliced


class CovarianceSearcher(PeakSearcher):

    def __init__(self, hwidth: int, FWHM: int, mode: Literal["uniform", "inverse", "normal"] = "inverse", label: str | None = None):
        self._hwidth = hwidth
        self._FWHM = FWHM
        self._mode = mode
        if label is None:
            label = f'CovarianceSearcher[F{self._FWHM}]'
        super().__init__(label=label)

    def __run__(self, spectra: list[Spectrum], *args, **kargs) -> Spectrum:
        searched = spectra[0].copy()
        covariance = self._covariance(searched)
        regions = self._split(covariance > 0)
        for region in regions:
            region.peaks.append({'location': round((region.left+region.right) / 2)})

        regions = self._filter_by_length(regions, searched)
        regions = self._filter_by_height(regions, searched)
        regions = self._expand(regions, searched)
        regions = self._merge(regions, searched)
        regions = self._filter_by_area(regions, searched)
        searched.regions = regions
        return searched

    def _covariance(self, spectrum: Spectrum) -> Spectrum:

        if self._hwidth > spectrum.counts.shape[0] - self._hwidth:
            raise ValueError(f"Half width {self._hwidth} is too large for spectrum with length {spectrum.counts.shape[0]}")

        covariances = np.zeros_like(spectrum)
        window = np.arange(-self._hwidth, self._hwidth+1)
        shape = np.exp(-4 * np.log(2) * (window/self._FWHM)**2)

        # pool = Pool(processes=4)
        # def func(i):
        #     windowed_spectrum = spectrum.counts[i-self._hwidth: i+self._hwidth+1]
        #     if self._mode == 'uniform':
        #         weight = np.ones(2*self._hwidth+1) / (2*self._hwidth+1)
        #     elif self._mode == 'inverse':
        #         weight = 1 / (windowed_spectrum + 1)
        #     elif self._mode == 'normal':  # with shape-counting weight
        #         weight = np.exp(-2*(windowed_spectrum/self._FWHM)**2) / windowed_spectrum
        #     else:
        #         raise ValueError(f"mode {self._mode} not supported.")

        #     variance = np.sum(weight)*np.sum(weight*shape**2) - np.sum(weight*shape)**2
        #     bias = np.power(np.sum(weight), 0.5) / variance**0.5
        #     covariance = np.sum(weight)*np.sum(weight*shape*windowed_spectrum) - np.sum(weight*windowed_spectrum)*np.sum(weight*shape)
        #     covariances[i] = covariance / variance / bias
        # pool.map(func, range(self._hwidth, spectrum.counts.shape[0]-self._hwidth))

        for i in range(self._hwidth, spectrum.counts.shape[0]-self._hwidth):
            windowed_spectrum = spectrum[i-self._hwidth: i+self._hwidth+1]
            if self._mode == 'uniform':
                weight = np.ones(2*self._hwidth+1)
            elif self._mode == 'inverse':
                weight = 1 / (windowed_spectrum + 1)
            elif self._mode == 'normal':  # with shape-counting weight
                
                weight = np.exp(-2*(windowed_spectrum.indexes/self._FWHM)**2) / windowed_spectrum
            else:
                raise ValueError(f"mode {self._mode} not supported.")

            variance = np.sum(weight)*np.sum(weight*shape**2) - np.sum(weight*shape)**2
            bias = np.power(np.sum(weight), 0.5) / variance**0.5
            covariance = np.sum(weight)*np.sum(weight*shape*windowed_spectrum) - np.sum(weight*windowed_spectrum)*np.sum(weight*shape)
            covariances[i] = covariance / variance / bias

        return Spectrum(covariances)


class SymmetricZeroAreaConvolutionSearcher(PeakSearcher):

    def __init__(self, hwidth: int, FWHM: int, func: Literal["gaussian", ] = "gaussian", label: str | None = None):
        self._hwidth = hwidth
        self._FWHM = FWHM
        self._func = func
        if label is None:
            label = f'SZACSearcher[F{self._FWHM}]'
        super().__init__(label=label)

    def __run__(self, spectra: list[Spectrum], *args, **kargs) -> Spectrum:
        searched = spectra[0].copy()
        convolved = self._convolve(searched)
        regions = self._split(convolved > 0)
        regions = self._filter_by_length(regions, searched)
        for region in regions:
            region.peaks.append({'location': round((region.left+region.right) / 2)})
        regions = self._filter_by_length(regions, searched)
        regions = self._expand(regions, searched)
        regions = self._merge(regions, searched)
        regions = self._filter_by_area(regions, searched)
        searched.regions = regions
        return searched

    def _convolve(self, spectrum: Spectrum) -> np.ndarray:
        indexes = np.arange(-self._hwidth, self._hwidth + 1)
        if self._func == "gaussian":
            weights = np.exp(-indexes**2/self._FWHM**2*4*np.log(2))
        elif self._func == "boxcar":
            weights = -1 * np.ones_like(indexes)
            weights[np.ceil(-(self._hwidth-1)/3): np.floor((self._hwidth-1)/3)] = 2
        elif self._func == "cauchy":
            weights = self._FWHM**2 / (self._FWHM**2 + 4*indexes**2)
        else:
            raise ValueError(f"func {self._func} not supported.")
        weights = weights - weights.mean()
        convolved = np.convolve(spectrum, weights, mode='same')
        return convolved


class SecondConvolutionPeakSearcher(PeakSearcher):

    def __init__(self, hwdith: int = 20, baseline: bool = False, label: str | None = None):
        self._hwidth = hwdith
        self._baseline = baseline
        self._indexes = np.arange(-self._hwidth, self._hwidth + 1)
        if label is None:
            label = f'SecondConvolutionSearcher[F{self._hwidth}]'
        super().__init__(label=label)

    def __run__(self, spectra: list[Spectrum], *args, **kargs) -> Spectrum:
        searched = spectra[0].copy()
        first_convol = self._firstconvol(spectra[0], std=self._hwidth/4)

        pos_regions = self._split(first_convol >= 0)
        neg_regions = self._split(first_convol < 0)
        pos_regions = self._filter_by_length(pos_regions, searched, 4)
        neg_regions = self._filter_by_length(neg_regions, searched, 4)

        regions = []
        for i, pos_region in enumerate(pos_regions):
            next_neg_regions = [neg_region for neg_region in neg_regions if neg_region.left >= pos_region.right]
            next_neg_region = next_neg_regions[0] if next_neg_regions else None
            if (next_neg_region is not None) and (i == len(pos_regions)-1 or next_neg_region.right < pos_regions[i+1].left):
                left = int((first_convol[pos_region.indexes]).argmax()) + pos_region.left
                right = int((first_convol[next_neg_region.indexes]).argmin()) + next_neg_region.left
                region = Region(left, right, peaks=[{'location': int((left+right) / 2)}])
                regions.append(region)

        regions = self._filter_by_length(regions, searched)
        regions = self._filter_by_height(regions, searched)
        regions = self._expand(regions, searched)

        merr, Merr = 0.5, 10
        stderrs = np.linspace(merr, Merr, 100)
        second_convols = [self._secondconvol(searched, stderr) for stderr in stderrs]
        for region in regions:
            best_stderr_index = np.array([second_convol[region.indexes].min() for second_convol in second_convols]).argmin()
            region.peaks[0]['stderr'] = stderrs[best_stderr_index]
            region.peaks[0]['center'] = int(second_convols[best_stderr_index][region.indexes].argmin()) + region.left
            region.peaks[0]['height'] = second_convols[best_stderr_index][region.indexes].min() * (- 3 * 3**0.5 / 2 / np.pi)
            region.peaks[0]['area'] = region.peaks[0]['stderr'] * region.peaks[0]['height'] * (2 * np.pi) ** 0.5

        regions = self._merge(regions, searched)
        regions = self._filter_by_area(regions, searched)
        searched.regions = regions

        if self._baseline:
            for region in searched.regions:
                centers = np.array([peak['center'] for peak in region.peaks])
                counts = np.array([searched[round(peak['center'])]-peak['height'] for peak in region.peaks])
                if len(centers) >= 2:
                    region.slope = ((centers * counts).sum() - centers.mean() * counts.sum()) / ((centers * centers).sum() - centers.mean() * centers.sum())
                    region.offset = counts.mean() - region.slope * centers.mean()
                else:
                    region.slope = 0
                    region.offset = searched[round(region.peaks[0]['center'])] - region.peaks[0]['height']
        else:
            for region in searched.regions:
                region.slope = 0
                region.offset = 0

        return searched

    @staticmethod
    def _gaussian_derivative(indexes, mean, std):
        return - np.exp(-(indexes - mean) ** 2 / (2 * std ** 2)) * (indexes - mean) / std ** 2

    def _firstconvol(self, spectrum: Spectrum, std: float) -> np.ndarray:
        weights = self._gaussian_derivative(self._indexes, 0, std)
        first_convol = np.convolve(spectrum.counts, weights, mode='same')
        return first_convol

    def _secondconvol(self, spectrum: Spectrum, stderr: float) -> np.ndarray:
        weights = self._gaussian_derivative(self._indexes, 0, stderr)
        first_convol = np.convolve(spectrum.counts, weights, mode='same')
        second_convol = np.convolve(first_convol, weights, mode='same')
        return second_convol


# if __name__ == '__main__':

#     import matplotlib.pyplot as plt
#     import Smoother
#     import Operator
#     import OtherOperator
#     import BasicOperator
#     from Spectrum import simuspecs

#     simu = simuspecs['doublepeak_slight']

#     cen = Smoother.CentroidSmoother(3)
#     strp = OtherOperator.SNIPStripper(10)
#     minus = BasicOperator.Stripper()

#     diff = DifferentialPeakSearcher(4, 3, 1, 2, 0.4)
#     gauss = GaussPeakSearcher(2, 0.05)
#     cov = CovarianceSearcher(2, 3)

#     smoothed = cen(cen(simu))
#     base = strp(cen(cen(simu)))
#     stripped = minus([cen(cen(simu)), base])

#     gaussspe = gauss(stripped)
#     diffspe = diff(stripped)
#     covspe = cov(stripped)

#     fig, axe = plt.subplots(3, 1)
#     for i, spe in enumerate([gaussspe, diffspe, covspe]):
#         spe.plot(axes=axe[i])
#         spe.plot_regions(axes=axe[i])
#         print(spe)
#     fig.savefig("compare.png")
