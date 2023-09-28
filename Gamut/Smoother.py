import numpy as np

from numpy.fft import fft, ifft
from pywt import wavedec, waverec

from Operator import Operator
from Spectrum import Spectrum


class CentroidSmoother(Operator):

    def __init__(self, order):
        self._order = order
        self._coefs = self._gencoefs()
        self._label = f'-CentroidSmoother[O{self._order}]'

    def __run__(self, spectrum):
        smthed = spectrum.copy()
        smthed.counts = np.convolve(spectrum.counts, self._coefs, 'same')
        smthed.label += self._label
        return smthed

    def _gencoefs(self):
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

    def __init__(self, order, hwidth):
        self._order = order
        self._hwidth = hwidth
        self._coefs = self._gencoefs()
        self._label = f'-SavitkySmoother[O{self._order}|HW{self._hwidth}]'

    def __run__(self, spectrum):
        smthed = spectrum.copy()
        smthed.counts = np.convolve(spectrum.counts, self._coefs, 'same')
        smthed.label += self._label
        return smthed

    def _gencoefs(self):
        mat_order, mat_width = np.meshgrid(np.arange(self._order+1),
                                           np.arange(-self._hwidth,
                                                     self._hwidth+1))
        A = mat_width ** mat_order
        M = np.matmul(np.linalg.inv(np.matmul(A.T, A)), A.T)
        return M[0, :]


class FourierSmoother(Operator):

    def __init__(self, mode, threshold):
        self._mode = mode
        self._threshold = threshold
        self._label = f'-FourierSmoother[{self._mode}|T{self._threshold:.2f}]'

    def __run__(self, spectrum):
        smthed = spectrum.copy()
        trsfed = fft(spectrum)
        trsfed = self._denoise(trsfed)
        smthed.counts = np.real(ifft(trsfed))
        smthed.label += self._label
        return smthed

    def _denoise(self, trsfed):
        if self._mode == 'low':
            dnsed = trsfed.copy()
            dnsed[int(self._threshold*trsfed.shape[0]):
                  int((1-self._threshold)*trsfed.shape[0])] = 0
        return dnsed


class WaveletSmoother(Operator):
    def __init__(self, wavelet, mode, order):
        self._wavelet = wavelet
        self._mode = mode
        self._order = order
        self._label = f'-WaveletSmoother[O{self._order}]'

    def __run__(self, spectrum):
        smthed = spectrum.copy()
        trsfed = wavedec(spectrum, wavelet=self._wavelet, level=self._order)
        trsfed = self._denoise(trsfed)
        smthed.counts = np.real(waverec(trsfed, wavelet=self._wavelet))
        smthed.label += self._label
        return smthed

    def _denoise(self, trsfed):
        sigma = np.median(np.abs(trsfed[-1]))
        dnsed = trsfed.copy()
        if self._mode == 'soft':
            thershold = sigma * np.sqrt(2 * np.log(len(trsfed)))
            for item in dnsed[-2:]:
                item[abs(item) < thershold] = 0
                item[item > thershold] -= thershold
                item[item < -thershold] += thershold
        return dnsed


class TranslationInvarianceWaveletSmoother(WaveletSmoother):

    def __init__(self, wavelet, mode, order, step=1):
        super().__init__(wavelet, mode, order)
        self._step = step
        self._label = f'-TIWaveletSmoother[O{self._order}]'

    def __run__(self, spectrum):
        smthed = spectrum.copy()
        self._wavelet_smoother = WaveletSmoother(self._wavelet, self._mode, self._order)
        trsfed_sum = np.zeros(spectrum.shape)
        times = spectrum.counts.shape[0] // self._step
        for i in range(times):
            trslted = Spectrum(np.roll(spectrum.counts, self._step*i))
            trsfed = self._wavelet_smoother(trslted)
            trsfed_sum += np.roll(trsfed.counts, -self._step*i)
        smthed.counts = trsfed_sum / times
        smthed.label += self._label
        return smthed


if __name__ == "__main__":

    import matplotlib.pyplot as plt
    from SimulatedSpectrum import SimulatedSpectrum

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