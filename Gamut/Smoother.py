import numpy as np

from typing import Literal
from numpy.fft import fft, ifft
from pywt import wavedec, waverec

from .Operator import Operator
from .Spectrum import Spectrum


class CentroidSmoother(Operator):

    def __init__(self, order: int = 2, label: str = None):
        '''
        The most basic smoother.

        :param order: smoothing order, also the half-width of smoothing window
        :label: label for the operator
        '''
        self._order = order
        self._coefs = self._gen_coefs()

        if label is None:
            label = f'CentroidSmoother[O{self._order}]'
        super().__init__(1, label)

    def __run__(self, spectra: Spectrum | list[Spectrum], *args, **kargs) -> Spectrum:
        smoothed = spectra[0].copy()
        padded = np.pad(spectra[0].data, (self._order, self._order), mode='reflect')
        smoothed[:] = np.convolve(padded, self._coefs, mode='same')[self._order: -self._order]
        return smoothed

    def _gen_coefs(self) -> np.ndarray:
        def coefs(order):
            if order == 0:
                return np.array([1])
            else:
                coef = coefs(order-1)
                coef1 = np.pad(coef, (2, 0), mode='constant',
                               constant_values=0)
                coef2 = np.pad(coef, (1, 1), mode='constant',
                               constant_values=0)
                coef3 = np.pad(coef, (0, 2), mode='constant',
                               constant_values=0)
            return coef1 / 4 + coef2 / 2 + coef3 / 4
        return coefs(self._order)


class SavitzkySmoother(Operator):

    def __init__(self, order: int = 2, hwidth: int = 3, label: str = None):
        '''
        The most widely used smoother, applied in commercial codes, such as GammaVision & Genie2000
        Reference:
            A. Savitzky and M. J. E. Golay, “Smoothing and Differentiation
            of Data by Simplified Least Squares Procedures.” Anal. Chem.,
            vol. 36, no. 8, pp. 1627–1639, Jul. 1964. [Online]. Available:
            http://dx.doi.org/10.1021/ac60214a047

        :param order: order of fitted polynomial
        :param hwidth: half width of smoothing window
        :param label: label for the operator
        '''
        self._order = order
        self._hwidth = hwidth
        self._coefs = self._gen_coefs()
        if label is None:
            label = f'SavitzkySmoother[O{self._order}]'
        super().__init__(1, label)

    def __run__(self, spectra: Spectrum | list[Spectrum], *args, **kargs) -> Spectrum:
        smoothed = spectra[0].copy()
        padded = np.pad(spectra[0].data, (self._hwidth, self._hwidth), mode='reflect')
        smoothed[:] = np.convolve(padded, self._coefs, mode='same')[self._hwidth: -self._hwidth]
        return smoothed

    def _gen_coefs(self) -> np.ndarray:
        mat_order, mat_width = np.meshgrid(np.arange(self._order+1),
                                           np.arange(-self._hwidth,
                                                     self._hwidth+1))
        mat_A = mat_width ** mat_order
        mat_M = np.matmul(np.linalg.inv(np.matmul(mat_A.T, mat_A)), mat_A.T)
        return mat_M[0, :]


class FourierSmoother(Operator):

    def __init__(self, mode: Literal["low", "high", "gauss"] = "low", threshold: float = 0.6, label: str = None):
        '''
        Smoother based on Fourier Transformation.
        Reference:


        :param order: order of fitted polynomial
        :param hwidth: half width of smoothing window
        :param label: label for the operator
        '''
        self._mode = mode
        self._threshold = threshold

        if label is None:
            label = f'FourierSmoother[{self._mode}]'
        super().__init__(1, label)

    def __run__(self, spectra: Spectrum | list[Spectrum], *args, **kargs) -> Spectrum:
        smoothed = spectra[0].copy()
        transformed = fft(spectra[0])
        transformed = self._denoise(transformed)
        smoothed[:] = np.real(ifft(transformed))
        return smoothed

    def _denoise(self, transformed: np.ndarray):
        if self._mode == 'low':
            denoised = transformed.copy()
            denoised[int(self._threshold*transformed.shape[0]):
                  int((1-self._threshold)*transformed.shape[0])] = 0
        return denoised


class WaveletSmoother(Operator):
    def __init__(self, wavelet: Literal[""], mode: Literal[""], order: int = 3):
        self._wavelet = wavelet
        self._mode = mode
        self._order = order
        if label is None:
            label = f'WaveletSmoother[O{self._order}]'
        super().__init__(1, label)

    def __run__(self, spectra: Spectrum | list[Spectrum], *args, **kargs) -> Spectrum:
        smoothed = spectra[0].copy()
        transformed = wavedec(spectra[0], wavelet=self._wavelet, level=self._order)
        transformed = self._denoise(transformed)
        smoothed[:] = np.real(waverec(transformed, wavelet=self._wavelet))
        return smoothed

    def _denoise(self, transformed: np.ndarray) -> np.ndarray:
        sigma = np.median(np.abs(transformed[-1]))
        denoised = transformed.copy()
        if self._mode == 'soft':
            thershold = sigma * np.sqrt(2 * np.log(transformed.shape[0]))
            for item in denoised[-2:]:
                item[abs(item) < thershold] = 0
                item[item > thershold] -= thershold
                item[item < -thershold] += thershold
        return denoised


class TranslationInvarianceWaveletSmoother(WaveletSmoother):

    def __init__(self, wavelet: Literal[""], mode: Literal[""], step: int = 1, order: int = 3):
        
        self._step = step
        if label is None:
            label = f'TIWaveletSmoother[O{self._order}]'
        super().__init__(wavelet, mode, order)

    def __run__(self, spectra: Spectrum | list[Spectrum], *args, **kargs) -> Spectrum:
        smoothed = spectra[0].copy()
        self._wavelet_smoother = WaveletSmoother(self._wavelet, self._mode, self._order)
        transformed_sum = np.zeros(spectrum.shape)
        times = spectra[0].shape[0] // self._step
        for i in range(times):
            translated = Spectrum(np.roll(spectrum, self._step * i))
            transformed = self._wavelet_smoother(translated)
            transformed_sum += np.roll(transformed, -self._step * i)
        smoothed[:] = transformed_sum / times
        smoothed.label += self._label
        return smoothed


if __name__ == "__main__":

    import matplotlib.pyplot as plt
    from Spectrum import SimulatedSpectrum

    spectrum = SimulatedSpectrum()

    centroid = CentroidSmoother(3)
    savitzky = SavitzkySmoother(3, 3)
    fourier = FourierSmoother('low', 0.2)
    wavelet = WaveletSmoother('db4', 'soft', 3)
    ti_wavelet = TranslationInvarianceWaveletSmoother('db4', 'soft', 3, 1)

    spe_cen = centroid(spectrum)
    spe_sav = savitzky(spectrum)
    spe_fou = fourier(spectrum)
    spe_wav = wavelet(spectrum)
    spe_ti_wav = ti_wavelet(spectrum)

    ax = spectrum.plot('.')
    ax_cen = spe_cen.plot()
    ax_sav = spe_sav.plot()
    ax_fou = spe_fou.plot()
    ax_wav = spe_wav.plot()
    ax_ti_wav = spe_ti_wav.plot()

    plt.ylim([10, 1000])
    plt.legend(fontsize=12)
    plt.show()

    # cen_smoothers = []
    # cen_spectrums = []
    # for i in range(1, 5):
    #     cen_smoothers.append(CentroidSmoother(i))
    #     cen_spectrums.append(cen_smoothers[-1](spectrum))
    #     cen_spectrums[-1].plot()
    # plt.ylim([10, 1000])
    # plt.legend(fontsize=6)
    # plt.show()