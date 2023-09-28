from Operator import Operator
from SimulatedSpectrum import SimulatedSpectrum
import numpy as np
from abc import abstractmethod
from utils import Differential, peak_fit

class PeakSearcher(Operator):

    @abstractmethod
    def __init__(self):
        pass

    def _split(self, data):
        slices = []
        start = None
        for i, d in enumerate(data):
            if d:
                if start is None:
                    start = i
            else:
                if start is not None:
                    slices.append((start, i-1))
                    start = None
        if start is not None:
            slices.append((start, i))
        return slices

    def _remove(self, peaks, threshold):
        return [(start, end) for start, end in peaks if end - start > threshold]

    def _correct(self, peaks, shift, spectrum):
        peaks_copy = []
        for start, end in peaks:
            start = max(start+shift, 0)
            end = min(end + shift, spectrum.shape[0] - 1)
            peaks_copy.append((start, end))
        return peaks_copy


class GaussPeakSearcher(PeakSearcher):
    '''
    if even o, o = o/2, g(i) = c(i+o) * c(i-o) / c(i+o+2) / c(i-o-2)
    if odd o, o = (o-1)/2 g(i) = c(i+o+1) * c(i-o) / c(i+o+3) / c(i-o-2)

    '''
    def __init__(self, order, threshold):
        self._order = order
        self._threshold = threshold
        self._label = f'-GaussPeakSearcher[O{self._order}]'

    def __run__(self, spectrum):
        searched = spectrum.copy()
        gauss = spectrum[2: -self._order] * spectrum[self._order: -2] \
            / spectrum[0: -self._order-2] / spectrum[self._order+2: ]
        ratio = (gauss - 1) * spectrum[2: -self._order]**0.5
        ratio = np.pad(ratio, (2, self._order), mode='edge')
        peaks = self._split(ratio > self._threshold)
        peaks = self._remove(peaks, 4)
        peaks = self._correct(peaks, int(self._order/2-1), spectrum)
        searched.attrbs['peaks'] = peaks
        searched.label += self._label
        return searched


class DifferentialPeakSearcher(PeakSearcher):

    def __init__(self, poly_order, hwidth, derive_order):
        self._poly_order = poly_order
        self._hwidth = hwidth
        self._derive_order = derive_order
        self._label = f'-DifferentialPeakSearcher[O{self._poly_order}]'

    def __run__(self, spectrum):
        searched = spectrum.copy()
        diff, _ = Differential(spectrum, self._poly_order, self._hwidth, self._derive_order)
        # return diff, diff[9:] * diff[:-9] < 0
        peaks = self._split(diff[9:] * diff[:-9] < 0)
        peaks = self._remove(peaks, 4)
        peaks = self._correct(peaks, 5, spectrum)
        searched.attrbs['peaks'] = peaks
        searched.label += self._label
        print(peaks)
        return searched


class PeakIdentifier(Operator):

    def __init__(self):
        pass

    def __run__(self, spectrum):
        identified = spectrum.copy()
        new_peaks = []
        for start, end in spectrum.attrbs['peaks']:
            windowed = spectrum[start: end+1]
            amplitude, centroid, stderr = peak_fit(windowed, 1)
            print(start, end, start+centroid, stderr)
            stderr = abs(stderr)
            if stderr < 6:
                new_peaks.append((max(int(start+centroid-3*stderr), 0), min(int(start+centroid+3*stderr), len(spectrum.counts)-1)))
        identified.attrbs['peaks'] = new_peaks
        print(new_peaks)
        return identified


class CovarianceSearcher(PeakSearcher):
    
    def __init__(self, hwidth, FWHM, mode):
        self._hwidth = hwidth
        self._FWHM = FWHM
        self._mode = mode
        self._label = f'-CovarianceSearcher[{self._mode}]'
    
    def __run__(self, spectrum):
        searched = spectrum.copy()
        cov = self._covariance(spectrum)
        peaks = self._split(cov > 0)
        peaks = self._remove(peaks, 2)
        searched.attrbs['peaks'] = peaks
        searched.label += self._label
        return searched

    def _covariance(self, spectrum):

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

        return covariances / variance / bias


if __name__ == '__main__':

    import matplotlib.pyplot as plt
    import Smoother
    import Operator
    import OtherOperator

    spe = SimulatedSpectrum()
    cen = Smoother.CentroidSmoother(2)
    strp = OtherOperator.SNIPstripper(5)
    diff = DifferentialPeakSearcher(3, 2, 1)
    gauss = GaussPeakSearcher(2, 0.1)
    cov = CovarianceSearcher(2, 3, 'gauss')
    ide = PeakIdentifier()

    spe = cen(spe)
    spe2 = spe.copy()
    spe = strp(spe)
    
    gaussspe = gauss(spe)
    diffspe = diff(spe)
    covspe = cov(spe)

    idegauss = ide(gaussspe)
    idediff = ide(diffspe)
    idecov = ide(covspe)
    idediff.counts = spe2.counts
    idegauss.counts = spe2.counts
    idecov.counts = spe2.counts


    fig, axe = plt.subplots(3, 1)
    # diffspe.plot(plot_peaks=True, ax=axe[0])
    # idespe.plot(plot_peaks=True, ax=axe[1])
    idegauss.plot(plot_peaks=True, ax=axe[0])
    idediff.plot(plot_peaks=True, ax=axe[1])
    idecov.plot(plot_peaks=True, ax=axe[2])
    plt.show()


    # gauss = GaussPeakSearcher(2, 0.1)
    # extend = BasicOperator.Extender(10)
    # slicer = BasicOperator.Slicer(10, -10)

    # pipe = Operator.OperatorPipe(cen, )
    # spe2 = pipe(spe)

    # spe3 = gauss(spe)

    # fig, axes = plt.subplots(3, 1, sharex=True)
    # spe.plot('.', ax=axes[0])
    # spe2.plot('.', ax=axes[1], plot_peaks=True)
    # ax = spe3.plot('.', ax=axes[2], plot_peaks=True)
    # # ax.twinx().plot(ratio, '.', label='ratio', color='red')
    # plt.show()

