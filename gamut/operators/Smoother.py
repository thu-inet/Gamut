import numpy as np
import scipy.sparse as sp

from typing import Literal
from numpy.fft import fft, ifft
from pywt import wavedec, waverec

from ..Operator import Operator
from ..Spectrum import Spectrum


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
        padded = np.pad(spectra[0].data, (self._order, self._order), mode='edge')  # coefficients of centroid method are symmetric, no reversion is required
        smoothed[:] = np.convolve(padded, self._coefs, mode='same')[self._order: -self._order]
        smoothed.propagate_covariance(self._coefs)
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

    def __init__(self, order: int = 2, hwidth: int = 2, label: str = None):
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
            label = f'SavitzkySmoother[O{self._order}W{2*self._hwidth+1}]'
        super().__init__(1, label)

    def __run__(self, spectra: Spectrum | list[Spectrum], *args, **kargs) -> Spectrum:
        smoothed = spectra[0].copy()
        padded = np.pad(spectra[0].data, (self._hwidth, self._hwidth), mode='reflect')
        smoothed[:] = np.convolve(padded, self._coefs, mode='same')[self._hwidth: -self._hwidth]
        smoothed.propagate_covariance(self._coefs)
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

    def _denoise(self, transformed: np.ndarray) -> np.ndarray:
        denoised = transformed.copy()
        if self._mode == 'low':
            denoised = transformed.copy()
            denoised[int(self._threshold*transformed.shape[0]):
                     int((1-self._threshold)*transformed.shape[0])] = 0
        return denoised


class WaveletSmoother(Operator):
    def __init__(self, wavelet: str | Literal["harr", "db1", "sym1", "coif1", "gaus1"] = 'sym3',
                 mode: Literal["hard", "soft", "quadratic-soft"] = "quadratic-soft",
                 threshold_mode: Literal["visushrink", "sqtwolog", "minimax"] = "visushrink",
                 order: int = 4, smooth_level: int=1, label: str | None = None):
        self._wavelet = wavelet
        self._mode = mode
        self._order = order
        self._threshold_mode = threshold_mode
        self._smooth_level = smooth_level
        if label is None:
            label = f'WaveletSmoother[O{self._order}M{self._mode[0]}T{self._threshold_mode[0]}]'
        super().__init__(1, label)

    def _shannon_entropy(self, transformed: list[np.ndarray]):
        SE_cA = - (transformed[0]**2 * np.log(transformed[0]**2)).sum()
        SE_cD = - (transformed[1]**2 * np.log(transformed[1]**2)).sum()
        return SE_cA, SE_cD

    def __run__(self, spectra: Spectrum | list[Spectrum], *args, **kargs) -> Spectrum:
        smoothed = spectra[0].copy()
        padded = np.pad(smoothed, pad_width=(100, 100), mode='edge')
        transformed = wavedec(padded, wavelet=self._wavelet, level=self._order, mode='smooth')
        transformed = self._denoise(transformed, padded)
        smoothed[:] = np.real(waverec(transformed, wavelet=self._wavelet))[100: -100]
        return smoothed

    def _denoise(self, transformed: list, smoothed: Spectrum) -> list:

        denoised = transformed.copy()
        sigma = np.median(np.abs(np.concatenate(transformed[1:])), axis=0) / 0.6745

        threshold = 0
        if self._threshold_mode == 'visushrink':
            threshold = sigma * (2 * np.log2(smoothed.shape[0])) ** 0.5
        elif self._threshold_mode == 'sqtwolog':
            threshold = (2 * np.log2(smoothed.shape[0])) ** 0.5
        elif self._threshold_mode == 'minimax':
            threshold = 0.3936 + 0.1829 * np.log2(smoothed.shape[0]) if smoothed.shape[0] > 32 else 0.0
            threshold = sigma * threshold
        elif self._threshold_mode == 'rigrsure':
            pass
        else:
            raise ValueError(f'Unknown threshold mode: {self._threshold_mode}')

        if self._mode == 'hard':
            for item in denoised[self._smooth_level:]:
                item[np.abs(item) < threshold] = 0
        elif self._mode == 'soft':
            for item in denoised[self._smooth_level:]:
                item[np.abs(item) < threshold] = 0
                indexes = np.abs(item) >= threshold
                item[indexes] = np.sign(item[indexes]) * np.maximum(np.abs(item[indexes]) - threshold, 0)
        elif self._mode == 'quadratic-soft':
            for item in denoised[self._smooth_level:]:
                item[abs(item) < threshold] = 0
                indexes = np.abs(item) >= threshold
                item[indexes] = np.sign(item[indexes]) * (item[indexes]**2 - threshold**2)**0.5
        else:
            raise ValueError(f'Unknown denoising mode: {self._mode}')
        return denoised


class TranslationInvarianceWaveletSmoother(WaveletSmoother):

    def __init__(self, wavelet: str | Literal["harr", "db2", "sym3", "gaus1", "dmey"] = 'sym3',
                 mode: Literal["hard", "soft", "quadratic-soft"] = 'quadratic-soft',
                 threshold_mode: Literal["visushrink", "sqtwolog", "minimax"] = "visushrink", step: int = 40, order: int = 3, label: str | None = None):
        self._step = step
        super().__init__(wavelet, mode, threshold_mode=threshold_mode, order=order, label=label)
        if label is None:
            self._label = f'TIWaveletSmoother[O{order}M{self._mode[0]}T{self._threshold_mode[0]}]'

    def __run__(self, spectra: Spectrum | list[Spectrum], *args, **kargs) -> Spectrum:
        smoothed = spectra[0].copy()
        transformed_sum = np.zeros(smoothed.shape)
        for i in range(-self._step, self._step+1):
            translated = Spectrum(np.roll(spectra[0], i))
            transformed = super().__run__([translated])
            transformed_sum += np.roll(transformed, -i)
        smoothed[15:-15] = transformed_sum[15:-15] / (2*self._step+1)
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
