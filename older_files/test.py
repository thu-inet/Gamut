import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import simu_spectrum

spectrum = simu_spectrum.simu_spectrum()
indexes = np.array([1, 2, 3])
counts = np.array([10**1.5, 40**1.5, 90**1.5]) / 1000

def Peak_fit(spectrum, region, npeak, mode='gauss'):
    indexes = np.arange(region[0], region[1]) 
    counts = spectrum[region[0]: region[1]]

    gauss = lambda x, h, x0, s: h * np.exp(-(x-x0)**2 / 2 / s**2)
    baseline = lambda x, k, b0: x * k + b0

    def fitfunc(x, *vars):
        y = 0
        for i in range(npeak):
            y += gauss(indexes, *vars[i*3: i*3+3])
        y += baseline(indexes, *vars[-2:])
        return y
    
    low = np.zeros(3*npeak+2)
    high = np.zeros(3*npeak+2)
    for i in range(npeak):
        low[i*3] = np.mean(sorted(counts)[:10])
        high[i*3] = np.mean(sorted(counts)[-10:])
        low[i*3+1] = 1
        high[i*3+1] = 10
        low[i*3+2] = region[0]
        high[i*3+2] = region[1]
    low[-2] = -10
    high[-2] = 10
    low[-1] = np.mean(sorted(counts)[:10])
    high[-1] = np.mean(sorted(counts)[-10:])
    vars, covs = curve_fit(fitfunc, indexes, counts, p0=(low+high)/2, bounds=(low, high))
    
    return vars
vars = Peak_fit(spectrum, [10, 40], 1)

# plt.plot(counts2)
# plt.plot(counts)
# plt.show()