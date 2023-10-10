from Operator import Operator
from Spectrum import Spectrum
import numpy as np
from utils import Differential, peak_fit
from PeakRegion import PeakRegion
from typing import Literal


class PeakSearcher(Operator):

    def __init__(self, label: str = None):
        super().__init__(1, label)

    def _split(self, data: np.ndarray) -> list[PeakRegion]:
        """
        Split a bool sequence into pieces of consecutive True.
        """
        peaks = []
        start = None
        for i, d in enumerate(data):
            if d:
                if start is None:
                    start = i
            else:
                if start is not None:
                    peaks.append(PeakRegion(start, i-1))
                    start = None
        if start is not None:
            peaks.append(PeakRegion(start, i-1))
        return peaks

    def _remove(self, peaks: list[PeakRegion], threshold: int):
        peaks_copy = []
        for peak in peaks:
            if peak.length >= threshold:
                peaks_copy.append(peak)
        return peaks_copy

    def _correct(self, peaks: list[PeakRegion], shift: int, spectrum: Spectrum):
        peaks_copy = []
        for peak in peaks:
            peak_copy = peak.copy()
            peak_copy.left = max(peak.left+shift, 0)
            peak_copy.right = min(peak.right + shift, spectrum.shape[0] - 1)
            peaks_copy.append(peak_copy)
        return peaks_copy


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
        peaks = self._split(padded > self._threshold)
        peaks = self._remove(peaks, 6)
        peaks = self._correct(peaks, round(self._order/2-1), searched)
        searched.peaks = peaks
        return searched

    def _transform(self, searched: Spectrum) -> Spectrum:
        padded = np.pad(searched, (2, self._order), mode='mean')
        gaussed = padded[2: -self._order] * padded[self._order: -2] / padded[0: -self._order-2] / padded[self._order+2: ]
        rectified = (gaussed - 1) * searched ** 0.5
        return Spectrum(rectified)


class DifferentialPeakSearcher(PeakSearcher):

    def __init__(self, poly_order: int, hwidth: int, derive_order: int, zero_width: int = 4, sigma: float = 0.7, label: str = None):
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
        zeros = np.zeros(differentiated.shape, dtype=bool)
        for i in range(self._zero_width, differentiated.shape[0]-self._zero_width):
            zeros[i] = np.all(differentiated[i-self._zero_width: i] > 0) and np.all(differentiated[i: i+self._zero_width] < 0)
        peaks = self._split(zeros)

        for peak in peaks:
            peak_hwidth = 1
            peak.left = max(0, peak.left - 2)
            peak.right = min(searched.shape[0] - 1, peak.right + 2)
            while True:
                peak.left = max(0, peak.left - 1)
                peak.right = min(searched.shape[0] - 1, peak.right + 1)
                total, baseline = searched.area_estimate(peak)
                if (total - baseline) / total  > self._sigma or ((total - baseline) / (total * peak_hwidth) ** 0.5  < 0.2):
                    break
                peak_hwidth += 1

        peaks = self._remove(peaks, 3)
        searched.peaks = peaks
        return searched

    def _transform(self, searched: Spectrum) -> np.ndarray:
        padded = np.pad(searched, (self._hwidth, self._hwidth), mode='edge')
        differentiated, _ = Differential(padded, self._poly_order, self._hwidth, self._derive_order)
        sliced = differentiated[self._hwidth: -self._hwidth]
        return Spectrum(sliced)


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
        peaks = self._split(covariance > 0)
        peaks = self._remove(peaks, 2)
        searched.peaks = peaks
        return searched

    def _covariance(self, spectrum: Spectrum) -> Spectrum:

        covariances = np.zeros(spectrum.counts.shape)
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
        spe.plot_peaks(axes=axe[i])
        print(spe)
    fig.savefig("compare.png")


