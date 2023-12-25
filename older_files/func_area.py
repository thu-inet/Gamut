import numpy as np
import matplotlib.pyplot as plt
import simu_spectrum
import func_smooth
from scipy.optimize import curve_fit

 
def Baseline_fit(spectrum, region, width=3, mode='linear'):
    index = np.concatenate((np.arange(region[0]-width, region[0]), np.arange(region[1], region[1]+width))) / region[1]
    counts = np.concatenate((spectrum[region[0]-width: region[0]], spectrum[region[1]: region[1]+width]))
    if mode == 'linear':
        basefunc = lambda x: np.vstack((x, np.ones_like(x))).T
    elif mode == 'quadratic':
        basefunc = lambda x: np.vstack((x**2, x, np.ones_like(x))).T
    elif mode == 'cubic':
        basefunc = lambda x: np.vstack((x**3, x**2, x, np.ones_like(x))).T
    elif mode == 'exponential':
        basefunc = lambda x: np.vstack((np.exp(x), np.ones_like(x))).T
    
    A = basefunc(index)
    b = np.linalg.inv(A.T @ A) @ A.T @ counts
    index_baseline = np.arange(region[0]-width, region[1]+width) 
    counts_baseline = basefunc(index_baseline / region[1] ) @ b 
    return index_baseline, counts_baseline

def Baseline_clip(spectrum, region):
    
    baseline = Baseline_fit(spectrum, region, width=3, mode='linear')[1]
    area_total = np.sum(spectrum[region[0]: region[1]])
    area_baseline = np.sum(baseline)

    return area_total - area_baseline


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
        low[i*3] = np.min(counts)
        high[i*3] = np.max(counts)
        low[i*3+1] = region[0]
        high[i*3+1] = region[1]
        low[i*3+2] = 1
        high[i*3+2] = 7
    low[-2] = -5
    high[-2] = 5
    low[-1] = np.mean(sorted(counts)[:10])
    high[-1] = np.mean(sorted(counts)[-10:])

    vars = (low+high)/2
    for i in range(npeak):
        vars[i*3+1] = (high[i*3+1]-low[i*3+1]) * (i+1)/(npeak+2) + low[i*3+1]

    vars, covs = curve_fit(fitfunc, indexes, counts, p0=vars, bounds=(low, high))

    return vars
    
if __name__ == '__main__':

    spectrum = simu_spectrum.simu_spectrum()
    # spectrum = func_smooth.Translation_Invariance_Wavelet(spectrum)
    # baseline = func_smooth.SNIP(spectrum, half_width=10)
    ind_left, ind_right = TPA(spectrum, [105, 115])
    # plt.plot(baseline)
    # plt.plot(spectrum)
    # for mode in ['linear', 'quadratic', 'cubic', 'exponential']:
    #     index_baseline, counts_baseline = Baseline_fit(baseline, [100, 120], width=7, mode=mode)
    #     plt.plot(index_baseline, counts_baseline, '--', label=mode)
    # plt.legend()
    # plt.show()

