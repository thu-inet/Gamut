import numpy as np


def gausspeak(centroid, stderror, area, channels):
    channels_peak = np.zeros(shape=channels.size)
    channels_peak[max(int(centroid-4*stderror), 0):min(int(centroid+4*stderror), channels.shape[0])] = \
        channels[max(int(centroid-4*stderror), 0):min(int(centroid+4*stderror), channels.shape[0])].copy()
    amplitude = area / stderror / (2*np.pi)**0.5
    return np.exp(- (channels - centroid) ** 2 / stderror ** 2 / 2) * amplitude


def baseline(channels, base_intensity=100, uncertainty=5, noise=0.5):
    base = -1 * \
        (channels / channels.shape[0] - 1)**2 * uncertainty + base_intensity
    noise = np.random.normal(loc=base, scale=base**noise, size=channels.shape)
    return noise


def simu_spectrum(peak_info=np.array([[30, 2.5, 4000], [70, 3, 2000], [110, 4, 1000],
                                      [150, 4.5, 400], [165, 4.3, 600]]), length=200):

    channels = np.arange(0, length)
    spectrum = baseline(channels)
    for i, item in enumerate(peak_info):
        peak = gausspeak(peak_info[i][0], peak_info[i]
                         [1], peak_info[i][2], channels)
        spectrum += peak
    return spectrum


if __name__ == '__main__':

    import matplotlib.pyplot as plt

    spectrum = simu_spectrum(length=1000)
    plt.plot(spectrum)
    plt.show()
