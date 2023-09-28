from Spectrum import Spectrum
import numpy as np


class SimulatedSpectrum(Spectrum):

    def __init__(self,
                 length=200,
                 peaks_info=[[30, 2.5, 4000], [70, 3, 2000], [110, 4, 1000],  [150, 4.5, 2400], [165, 4.3, 1600]],
                 base_intensity=100,
                 base_amplitude=100,
                 base_function=lambda x: x,
                 label='simulated'):
        self.length = length
        self.peaks_info = peaks_info
        self.base_intensity = base_intensity
        self.base_amplitude = base_amplitude
        self.base_function = base_function

        peaks_index = []
        for peak_info in peaks_info:
            peaks_index.append((int(peak_info[0]-peak_info[1]*3),
                                int(peak_info[0]+peak_info[1]*3+1)))
        self.base, self.peak = self._genspec()
        super().__init__(self.base+self.peak, label, peaks=peaks_index)

    def _genspec(self):
        base = self.base_function(np.linspace(0, 1, self.length))
        base = (base - base.min()) / base.ptp()
        base = base * self.base_amplitude + self.base_intensity
        base = np.random.normal(loc=base, scale=base**0.5, size=base.shape)

        peak = np.zeros(shape=base.shape)
        for peak_info in self.peaks_info:
            peak += self._genpeak(*peak_info)
        return base, peak

    def _genpeak(self, centroid, stderror, area):
        channels = np.arange(0, self.length)
        amplitude = area / stderror / (2*np.pi)**0.5
        return amplitude * np.exp(- (channels-centroid)**2 / stderror**2 / 2)

if __name__ == '__main__':

    import matplotlib.pyplot as plt

    spe = SimulatedSpectrum()
    spe.plot()
    plt.show()