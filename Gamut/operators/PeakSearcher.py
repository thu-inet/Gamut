import numpy as np

from typing import Literal

from ..utils import Differential
from ..PeakRegion import Region

from ..Operator import Operator
from ..Spectrum import Spectrum


class PeakSearcher(Operator):

    def __init__(self, label: str = None):
        super().__init__(1, label)

    def _split(self, data: np.ndarray) -> list[Region]:
        """
        Split a bool sequence into pieces of consecutive True.
        """
        regions = []
        start = None
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

    def _correct(self, regions: list[Region], shift: int, spectrum: Spectrum):
        regions_copy = []
        for peak in regions:
            peak_copy = peak.copy()
            peak_copy.left = max(peak.left+shift, 0)
            peak_copy.right = min(peak.right + shift, spectrum.shape[0] - 1)
            regions_copy.append(peak_copy)
        return regions_copy

    def _remove_length(self, regions: list[Region], length: int):
        regions_copy = []
        for region in regions:
            if region.length >= length:
                regions_copy.append(region)
        return regions_copy

    def _remove_height(self, regions: list[Region], height: int, spectrum: Spectrum):
        regions_copy = []
        for region in regions:
            if spectrum[region.indexes].max()-spectrum[region.indexes].mean() >= 0.582 * height:
                regions_copy.append(region)
        return regions_copy

    def _expand(self, regions: list[Region], spectrum: Spectrum):
        for region in regions:
            while True:
                total, baseline = spectrum.area_estimate(region)
                area = total - baseline
                expanded_region = [Region(max(region.left-1, 0), min(region.right+0, spectrum.length-1)),
                                   Region(max(region.left-0, 0), min(region.right+1, spectrum.length-1)),
                                   Region(max(region.left-1, 0), min(region.right+1, spectrum.length-1)),
                                   Region(max(region.left-2, 0), min(region.right+1, spectrum.length-1)),
                                   Region(max(region.left-1, 0), min(region.right+2, spectrum.length-1)),
                                   Region(max(region.left-2, 0), min(region.right+2, spectrum.length-1))]
                expanded_area = [spectrum.area_estimate(region2)[0]-spectrum.area_estimate(region2)[1] for region2 in expanded_region]
                expanded_region = np.array(expanded_region)
                if any(expanded_area > area):
                    region2 = expanded_region[np.argmax(expanded_area)]
                    region.left = region2.left
                    region.right = region2.right
                else:
                    break
        return regions

    def _remove_area(self, regions: list[Region], area_ratio: float, spectrum: Spectrum):
        regions_copy = []
        for region in regions:
            total, baseline = spectrum.area_estimate(region)
            if (total - baseline) / total >= area_ratio:
                regions_copy.append(region)
        return regions_copy

    def _merge(self, regions: list[Region]):
        max_index = max([region.right for region in regions])
        united_indexes = np.zeros(max_index+1)
        united_peaks = []
        for region in regions:
            united_indexes[region.indexes] = 1
            united_peaks.extend(region.peaks)
        merged_regions = self._split(united_indexes)
        for region in merged_regions:
            region.peaks = [peak for peak in united_peaks if (peak['location'] <= region.right) and (peak['location'] >= region.left)]
        return merged_regions


class GaussPeakSearcher(PeakSearcher):
    '''
    if even o, o = o/2, g(i) = c(i+o) * c(i-o) / c(i+o+2) / c(i-o-2)
    if odd o, o = (o-1)/2 g(i) = c(i+o+1) * c(i-o) / c(i+o+3) / c(i-o-2)

    '''
    def __init__(self, order: int, threshold: float, label: str = None):
        self._order = order
        self._threshold = threshold
        if label is None:
            label = f'GaussPeakSearcher[O{self._order}]'
        super().__init__(label)

    def __run__(self, spectra: Spectrum | list[Spectrum], *args, **kargs) -> Spectrum:
        searched = spectra[0].copy()
        padded = self._transform(searched)
        regions = self._split(padded > self._threshold)
        regions = self._correct(regions, round(self._order/2-1), searched)
        for region in regions:
            region.peaks.append({'location': int((region.left+region.right) / 2)})
        regions = self._remove_length(regions, 4)
        regions = self._expand(regions, searched)
        regions = self._merge(searched.regions + regions)
        searched.regions = regions
        return searched

    def _transform(self, searched: Spectrum) -> np.ndarray:
        padded = np.maximum(np.pad(searched, (2, self._order), mode='constant'), 1)
        gaussed = padded[2: -self._order] / padded[0: -self._order-2] / padded[self._order+2:] * padded[self._order: -2]
        rectified = (gaussed - 1) * searched ** 0.5
        return rectified


class DifferentialPeakSearcher(PeakSearcher):

    def __init__(self, poly_order: int = 1, hwidth: int = 1, derive_order: int = 1, zero_width: int = 4, sigma: float = 0.1, label: str = None):
        self._poly_order = poly_order
        self._hwidth = hwidth
        self._derive_order = derive_order
        self._zero_width = zero_width
        self._sigma = sigma
        if label is None:
            label = f'DifferentialPeakSearcher[O{self._poly_order}]'
        super().__init__(label)

    def __run__(self, spectra: Spectrum | list[Spectrum], *args, **kargs) -> Spectrum:
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
        regions = self._remove_length(regions, 4)
        regions = self._remove_height(regions, 10, searched)
        regions = self._remove_area(regions, self._sigma, searched)
        regions = self._expand(regions, searched)
        regions = self._merge(searched.regions + regions)
        searched.regions = regions
        return searched

    def _transform(self, searched: Spectrum) -> np.ndarray:
        padded = np.pad(searched, (self._hwidth, self._hwidth), mode='edge')
        differentiated, _ = Differential(padded, self._poly_order, self._hwidth, self._derive_order)
        sliced = differentiated[self._hwidth: -self._hwidth]
        return sliced


class CovarianceSearcher(PeakSearcher):

    def __init__(self, hwidth: int, FWHM: int, mode: Literal["uniform", "inverse", "normal"] = "inverse", label: str = None):
        self._hwidth = hwidth
        self._FWHM = FWHM
        self._mode = mode
        if label is None:
            label = f'CovarianceSearcher[F{self._FWHM}]'
        super().__init__(label)

    def __run__(self, spectra: Spectrum | list[Spectrum], *args, **kargs) -> Spectrum:
        searched = spectra[0].copy()
        covariance = self._covariance(searched)
        regions = self._split(covariance > 0)
        for region in regions:
            region.peaks.append({'location': round((region.left+region.right) / 2)})
        regions = self._remove_length(regions, 4)
        regions = self._expand(regions, searched)
        regions = self._merge(searched.regions + regions)
        searched.regions = regions
        return searched

    def _covariance(self, spectrum: Spectrum) -> Spectrum:

        # padded_spectrum = np.pad(spectrum, (self._hwidth, self._hwidth), mode='edge')
        covariances = np.zeros_like(spectrum)
        for i in range(self._hwidth, spectrum.counts.shape[0]-self._hwidth):
            windowed_spectrum = spectrum.counts[i-self._hwidth: i+self._hwidth+1]
            if self._mode == 'uniform':
                weight = np.ones(2*self._hwidth+1) / (2*self._hwidth+1)
            elif self._mode == 'inverse':
                weight = 1 / windowed_spectrum
            elif self._mode == 'normal':  # with shape-counting weight
                weight = np.exp(-2*(windowed_spectrum/self._FWHM)**4) / windowed_spectrum
            else:
                weight = np.ones(windowed_spectrum.shape)
            window = np.arange(-self._hwidth, self._hwidth+1)
            shape = np.exp(-4 * np.log(2) * (window/self._FWHM)**2)
            variance = sum(weight)*sum(weight*shape**2)-sum(weight*shape)**2
            bias = sum(weight)**0.5 / variance**0.5
            covariance = sum(weight)*sum(weight*shape*windowed_spectrum)-sum(weight*windowed_spectrum)*sum(weight*shape)
            covariances[i] = covariance
        return Spectrum(covariances / variance / bias)


class SymmetricZeroAreaConvolutionSearcher(PeakSearcher):

    def __init__(self, hwidth: int, FWHM: int, func: Literal["gaussian", ] = "gaussian", label: str = None):
        self._hwidth = hwidth
        self._FWHM = FWHM
        self._func = func
        if label is None:
            label = f'SZACSearcher[F{self._FWHM}]'
        super().__init__(label)

    def __run__(self, spectra: Spectrum | list[Spectrum], *args, **kargs) -> Spectrum:
        searched = spectra[0].copy()
        convolved = self._convolve(searched)
        regions = self._split(convolved > 0)
        regions = self._remove_length(regions, 4)
        for region in regions:
            region.peaks.append({'location': round((region.left+region.right) / 2)})
        regions = self._remove_length(regions, 4)
        regions = self._expand(regions, searched)
        regions = self._merge(searched.regions + regions)
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
        weights = weights - weights.mean()
        convolved = np.convolve(spectrum, weights, mode='same')
        return convolved


class SecondConvolutionPeakSearcher(PeakSearcher):

    def __init__(self, hwdith: int = 20, sigma: float = 0.1, baseline: bool = False, label: str = None):
        self._hwidth = hwdith
        self._sigma = sigma
        self._baseline = baseline
        self._indexes = np.arange(-self._hwidth, self._hwidth + 1)
        if label is None:
            label = f'SecondConvolutionSearcher[F{self._hwidth}]'
        super().__init__(label)

    def __run__(self, spectra: Spectrum | list[Spectrum], *args, **kargs) -> Spectrum:
        searched = spectra[0].copy()
        first_convol = self._firstconvol(spectra[0], std=self._hwidth/4)

        pos_regions = self._split(first_convol >= 0)
        neg_regions = self._split(first_convol <= 0)
        pos_regions = self._remove_length(pos_regions, 2)
        neg_regions = self._remove_length(neg_regions, 2)

        regions = []
        for i, pos_region in enumerate(pos_regions):
            next_neg_region = [neg_region for neg_region in neg_regions if neg_region.left >= pos_region.right][0]
            if (next_neg_region is not None) and (i == len(pos_regions)-1 or next_neg_region.right < pos_regions[i+1].left):
                left = (first_convol[pos_region.indexes]).argmax() + pos_region.left
                right = (first_convol[next_neg_region.indexes]).argmin() + next_neg_region.left
                region = Region(left, right, peaks=[{'location': int((left+right) / 2)}])
                regions.append(region)
        regions = self._remove_length(regions, 4)
        regions = self._remove_height(regions, 10, searched)
        regions = self._remove_area(regions, self._sigma, searched)
        regions = self._expand(regions, searched)

        merr, Merr = 1, 10
        stderrs = np.linspace(merr, Merr, 300)
        second_convols = [self._secondconvol(searched, stderr) for stderr in stderrs]
        for region in regions:
            # splits = region.split_peaks()
            best_stderr_index = np.array([second_convol[region.indexes].min() for second_convol in second_convols]).argmin()
            region.peaks[0]['stderr'] = stderrs[best_stderr_index]
            region.peaks[0]['center'] = second_convols[best_stderr_index][region.indexes].argmin() + region.left
            region.peaks[0]['height'] = second_convols[best_stderr_index][region.indexes].min() * (- 3 * 3**0.5 / 2 / np.pi)
            region.peaks[0]['area'] = region.peaks[0]['stderr'] * region.peaks[0]['height'] * (2 * np.pi) ** 0.5

        regions = self._merge(regions)
        searched.regions = regions

        if self._baseline:
            for region in searched.regions:
                centers = np.array([peak['center'] for peak in region.peaks])
                counts = np.array([searched[round(peak['center'])] for peak in region.peaks])
                if len(centers) >= 2:
                    region.slope = ((centers * counts).sum() - centers.mean() * counts.sum()) / (centers * centers).sum() - centers.mean() * centers.sum()
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

    def _firstconvol(self, spectrum: Spectrum, std: int) -> np.ndarray:
        weights = self._gaussian_derivative(self._indexes, 0, std)
        first_convol = np.convolve(spectrum.counts, weights, mode='same')
        return first_convol

    def _secondconvol(self, spectrum: Spectrum, stderr: int) -> np.ndarray:
        weights = self._gaussian_derivative(self._indexes, 0, stderr)
        first_convol = np.convolve(spectrum.counts, weights, mode='same')
        second_convol = np.convolve(first_convol, weights, mode='same')
        return second_convol


if __name__ == '__main__':

    import matplotlib.pyplot as plt
    import Smoother
    import Operator
    import OtherOperator
    import BasicOperator
    from Spectrum import simuspecs

    simu = simuspecs['doublepeak_slight']

    cen = Smoother.CentroidSmoother(3)
    strp = OtherOperator.SNIPStripper(10)
    minus = BasicOperator.Stripper()

    diff = DifferentialPeakSearcher(4, 3, 1, 2, 0.4)
    gauss = GaussPeakSearcher(2, 0.05)
    cov = CovarianceSearcher(2, 3)

    smoothed = cen(cen(simu))
    base = strp(cen(cen(simu)))
    stripped = minus([cen(cen(simu)), base])

    gaussspe = gauss(stripped)
    diffspe = diff(stripped)
    covspe = cov(stripped)

    fig, axe = plt.subplots(3, 1)
    for i, spe in enumerate([gaussspe, diffspe, covspe]):
        spe.plot(axes=axe[i])
        spe.plot_regions(axes=axe[i])
        print(spe)
    fig.savefig("compare.png")
