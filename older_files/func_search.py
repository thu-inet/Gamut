import class_objects
import numpy as np
from scipy.signal import find_peaks, argrelextrema

def list_split(l):
    '''
    Auxiliary function to split index lists into consective groups.
    Example: [1,2,3,6,7,8,9,12,13,14] -> [1,2,3] [6,7,8,9] [12,13,14]

    :param l: list or 1-d array to split
    :return splits: 2d array of splitted list
    '''
    # l = np.array(l, dtype=np.int16)
    splits = []
    index_start = 0
    index_end = 0
    
    while True:
        incre = l[min(index_end+1, len(l)-1)] - l[index_end]
        if  incre == 1:
            index_end += 1
        elif incre > 1:  
            splits.append(l[index_start: index_end+1])
            index_start = index_end + 1
            index_end = index_start
        elif incre == 0:
            splits.append(l[index_start: index_end+1])
            return splits
        else:
            raise TypeError("Not increasing list")

def threshold_split(spectrum, threshold):
    '''
    Wrapper of list_split, only return indexes
    Example: [1,2,3,6,7,1,2,3,6,7,8], threshold=4 -> [[3,4], [8,9,10]]

    :param spectrum: spectrum counts
    :param threshold: threshold to split

    :return splits: 2d array of splitted list
    '''
    if not isinstance(spectrum, np.ndarray):
        spectrum = np.array(spectrum)
    return list_split(np.where(spectrum >= threshold)[0])
    
def locate_peak(counts, peak_region, mode="median"):
    '''
    Locate the peak position in the peak region. Three modes are available: 
        -median: median of peak region
        -max: channel of maximum counts of peak region
        -weight: weighted average of peak region        
    
    :param counts: spectrum counts
    :param peak_region: peak region, (left, right)
    :param mode["median"]: mode of peak locating
    :return index: index of peak position
    '''
    windowed_counts = counts[peak_region[0]: peak_region[-1]]
    if mode == "max":
        index = windowed_counts.index(windowed_counts.max())
    elif mode == 'median':
        index = int((peak_region[0]+peak_region[-1]) * 0.5)
    elif mode == "weight":
        index = int(sum(windowed_counts*np.arange(windowed_counts.shape[0]))) + peak_region[0]
    
    return index

def locate_side_extrema(counts, peak):
    '''
    Locate the left and right side extrema of the given peak position.

    :param counts: spectrum counts
    :param peak: peak position

    :return (left, right): left and right side extrema    
    '''
    try:
        left = argrelextrema(counts[:peak], np.less)[0][-1]
    except:
        left = peak
    try:
        right =argrelextrema(counts[peak:], np.less)[0][0] + peak + 1
    except:
        right = peak + 1
    return (left, right)

def transform(coefs, counts):
    '''
    Wrapper of np.convolve, return the convolution of counts and coefs.
    '''
    return np.convolve(coefs[::-1], counts, mode='same')

def SNIP(counts, order=15, mode='positive'):
    '''
    Basic SNIP algorithm to clip the baseline component.
    (1) sqrt-log-log transform to compress variation range of spectrum counts
    conv =  log( log( (counts+1)^0.5+1 ) + 1 )
    (2) iteration to flatten the spectrum with symmetric points
    conv[i] = min(conv[i], 0.5*(conv[i-m]+conv[i+m]))
    (3) inverse transform to reestablish the spectrum
    baseline = exp( exp(conv)-1 )**2 - 1
    
    :param counts: spectrum counts
    :param order[15]: time to flatten the spectrum, and maximum distance from center point to alternative points
    :param mode["positive"]: iteration direction, "positive" (m=1->order) or "negative" (m=order->1)
    
    :return baseline: clipped baseline
    '''
    conv = np.log( np.log( (counts+1)**0.5+1 ) + 1 )
    for mi in range(1, order+1):
        mi = order+1 - mi if mode == 'positive' else mi
        conv[mi: -mi] = np.minimum(conv[mi: -mi], (conv[0: -2*mi]+conv[2*mi: ])/2)
    baseline = np.exp( np.exp(conv)-1 )**2 - 1
    
    return baseline

def plot_peaks(peak_list, counts, ax=None):
    
    if ax == None:
        ax = plt.gca()
        
    for i, peak in enumerate(peak_list):
        left, right = peak
        ax.plot(np.arange(left, right+1), counts[left: right+1], color='steelblue')
        
    return ax

def Hilbert(counts, order):
    '''
    Hilbert peak-searching algorithm:
    
    '''
    coefs_hilbert = np.arange(-order, order+1)
    coefs_hilbert[order] = 1
    coefs_hilbert = 1 / coefs_hilbert
    coefs_hilbert[order] = 0
    hilbert = np.convolve(coefs_hilbert[::-1], counts, mode='same')
    return hilbert, coefs_hilbert
    
def Gauss(counts, order=2):
    
    '''
    Gauss multiplication factor peak-searching algorithm:
        G[i] = c[i]*c[i+o-2] / c[i-2]*c[i+o], where o is the order
    Assuming G[i] ~ Aexp(-(i-i0)^2)+B near the peak region and o=2
        G[i0] = (A+B)^2 / (Aexp(-4)+B)^2 
        G should be bigger than 1 for peaks and equal to 1 for baseline
    In typycal spectrum, sigma ~ 3, exp ~ 0.8
        G[i0] = 1.22(A=B) 1.38(A=4B) 1.44(A=10B)
    The prominence is defined as the multiplication of G and noise 
    to tolerate possible fluctuation.
        P[i] = (G[i]-1) * (A+B)**0.5
    The G>1 region should be consecutive, it should be splited into pieces.
    and only pieces wider than 3 are considered as a peak
    
    :param counts: spectrum counts
    :param order[2]: order of Gauss peak-searching algorithm
    
    :return prominence_spectrum
    '''
    
    counts += 1
    gauss_spectrum = counts[2: -order] * counts[order: -2]  / counts[0: -order-2] / counts[2+order: ]
            
    prominence_spectrum = ( gauss_spectrum - 1 ) * counts[2: -order]**0.5 
    prominence_spectrum = np.pad(prominence_spectrum, pad_width=(2, order), mode='edge')

    return prominence_spectrum

def Differential(counts, order=2, half_width=3, derive_order=1):
    '''
    Numerical differential peak-searching method based on Savitzy-Colay fitting method.
    Equation: 
        Ab + E = y
        -b = [b0, b1, ..., bi,..., bn], fitted polynomial
        -A = [[1   m     m^2    ,...,       m^n] 
              [1  m-1  (m-1)^2  ,...,   (m-1)^n]
              [1   i     i^2    ,...,       i^n]
              [1  -m   (-m)^2   ,...,    (-m)^n]], index for fitted polynomial
        -y = [y(m), y(m-1),..., y(0),..., y(-m)], windowed spectrum
        -m = half_width of transformation window
        -n = order of fitted polynomial
        -bi = polynomial coefficients
    Solution: 
        -b = (A.t*A).I * A.t * y
        -y(fitted) = Ab = A * (A.t*A).I *A.t * y = A * M * y
        -y(0)(kth derivative)(fitted) 
            = y(0)(fitted)(kth derivative)
            = ( b0 + b1 * m + ... + bn * m^n )(kth derivative)|(m=0)
            = k! * bk + (k+1)!/1! * b(k+1) * m + ... + n!/(n-k)! * bn * m^(n-k)|(m=0)
            = k! * bk
            = K! * (M * y)[kth element] 
            = k! * (M * y)[k] 
        * M = (A.t*A).I *A.t, b = M * y
    
    :param counts: spectrum counts
    :param order[2]: polymonial order of target fitted function
    :param half_width[3]: halfwidth of transform window
    :param derive_order[1]: derivative order of numerical differential
    
    :return diff: target numerical differential
    :return coefs: coefficient vector of transform 
    '''
    mat_order, mat_width = np.meshgrid(np.arange(order+1), np.arange(-half_width, half_width+1))
    A = mat_width ** mat_order
    M = np.matmul( np.linalg.inv(np.matmul(A.T, A)), A.T )
    diff = np.zeros(counts.shape)
    if derive_order == 0:
        for i in range(half_width, counts.shape[0]-half_width):
            diff[i] = np.matmul( M, counts[i-half_width: i+half_width+1] )[derive_order]
        coefs = M[derive_order]
    else:
        for i in range(half_width, counts.shape[0]-half_width):
            diff[i] = np.matmul( M, counts[i-half_width: i+half_width+1] )[derive_order] * derive_order
        coefs = M[derive_order] * derive_order
    return diff, coefs
    

def Covariance(counts, half_width=8, FWHM=8, weight_mode='uniform'):
    
    '''
    Covariance peak-searching method, proposed by H.P.Blok in 1975, see ref: https://doi.org/10.1016/0029-554X(75)90523-6
    (1) Least square fit of spectrum counts to shape function
    yi+j = Ai * Sj + bi
    Sj is shape function, Ai and bi are fitted parameter at channel i.
    (2) Calculate the optimal paramters
    bi  = - sum[j](Ai * Sj) / N, j=-N~N
    Ai  = sum[j](Si * yi+j) - Si_avg * sum[j](yi+j) / sum[j](Si**2) - (sum[j](Si))**2
        = sum[j]((Si-S_avg) * yi+j) / sum[j](Si**2) - (sum[j](Si))**2
        = sum[j]((Si-S_avg) * (yi+j-yi+j_avg)) / sum[j](Si-Si_avg)**2
        = convariance between Si and yi / variance of Si
    (3) Add weights to Si and yi 
    Ai  = sum[j](Si * yi+j) - Si_avg * sum[j](yi+j) / sum[j](Si**2) - (sum[j](Si))**2
        = sum[j](gi) * sum[j](gi * Si * yi+j) - sum[i](gi * Si) * sum[j](gi * yi+j) 
        / sum[j](gi) * sum[j](gi * Si**2) - (sum[j](gi * Si))**2
    dAi =sum[j](gi)**0.5 / variance**0.5
    (4) Calculate final prominence
    pi = Ai / dAi = covariance / variance / bias
    
    :param counts: spectrum counts
    :param half_width[8]: halfwidth of convolution window
    :param FWHM[8]: estimated FWHM    
    :param weight_mode["uniform"]:  weight mode, "uniform", "inverse": gi=1/yi+j, "normal": gi=exp(-i**4/FWHM**4)/yi+j
    
    :return prominence
    '''
    
    if weight_mode == 'uniform': 
        weight = np.ones(2*half_width+1) / (2*half_width+1)
    elif weight_mode == 'inverse': 
        weight = 1 / windowed_spectrum
    elif weight_mode == 'normal': # with shape-counting weight
        weight = np.exp(-2*(window/FWHM)**4) / windowed_spectrum
    else:
        weight = np.ones(windowed_spectrum.shape)
    
    window = np.arange(-half_width, half_width+1)
    shape = np.exp( -4 * np.log(2) * (window/FWHM)**2 )
    variance = sum(weight)*sum(weight*shape**2)-sum(weight*shape)**2 
    bias = sum(weight)**0.5 / variance**0.5
    
    covariances = np.zeros(counts.shape)  
    for i in range(half_width, counts.shape[0]-half_width):
        windowed_spectrum = counts[i-half_width: i+half_width+1]
        covariance = sum(weight)*sum(weight*shape*windowed_spectrum)-sum(weight*windowed_spectrum)*sum(weight*shape) 
        covariances[i] = covariance 
    
    return covariances / variance / bias 

def SZAC(counts, half_width=10, scale=4, mode="gauss"):
    
    '''
    Sysmetric Zero Area Convolution peak-searching method \n
    (1) Convolution function so that conv1 = 0 for baseline, conv1!=0 for peak region 
    sum[j](Sj)=0 \n
    (2) Convolute the spectrum with kernel function \n
    conv1i = sum[j](yi+j * Sj) \n
    conv2i = sum[j](yi+j * Si**2) \n
    (3) Calculate the prominence \n
    pi = conv1i / conv2i \n
    
    :param counts: spectrum counts 
    :param half_width[5]: halfwidth of convolution windows 
    :param scale[4]: characteristic scale of convolution kernel
    :param mode["gauss"]: shape-like functions, "gauss", "rectangle", "cauchy"

    :return prominence     
    '''
    if mode == 'rectangle':
        func = lambda x: 1 if (x<1) and (x>=-1) else 0
    elif mode == 'gauss':
        func = lambda x: np.exp( -x**2 )  
    elif func == 'cauchy':
        func = lambda x: 1 / (1+x**2)
    
    window = np.arange(-half_width, half_width+1)
    coefs = func( window/scale ) 
    coefs = coefs - coefs.mean()
    trans1 = transform(coefs, counts)
    trans2 = transform(coefs**2, counts)
    
    return trans1 / trans2



if __name__ == '__main__':

    import matplotlib.pyplot as plt
    import simu_spectrum
    counts = simu_spectrum.simu_spectrum(length=300)
    
    diff, _ = Differential(counts)
    diff, _ = Differential(counts)

    # cov = Covariance(counts, 10, 15)
    # cov = SZAC(counts, 8, 3, 'Gauss')
    plt.plot(diff)
    plt.twinx()
    plt.plot(counts, linestyle='-', color='black')
    # plot_peaks(categorize(cov, 1.2, 1), counts)
    # for peak in list_peaks:
    #     centroid, left, right = peak.centroid, peak.left, peak.right
        # plt.plot(np.arange(left, min(right+1, counts.shape[0])), counts[left: right+1], 'x')
        # plt.vlines([left, right], ymin=0, ymax=1500)
    plt.show()
    