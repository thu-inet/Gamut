import numpy as np
import pywt

def centroid_coefficients(order: int):
    '''
    Calculate centroid smoothing coefficients by recursion.
    
    :param i: order of centroid method
    :return : coefficient array
    '''
    if order == 0:
        return np.array([1])
    else:
        coef = centroid_coefficients(order-1)
        coef1 = np.pad(coef, pad_width=(2, 0), mode='constant', constant_values=0)
        coef2 = np.pad(coef, pad_width=(1, 1), mode='constant', constant_values=0)
        coef3 = np.pad(coef, pad_width=(0, 2), mode='constant', constant_values=0)
        return coef1/4 + coef2/2 + coef3/4
    

def Centroid(spectrum: np.array, half_width: int = 2):
    '''
    centroid smoothing method
    
    :param spectrum: analysing spectrum
    :param half_width: halfwidth of transform window
    :return smoothed: smoothed spectrum
    '''
    N = spectrum.shape[0]
    coefs = centroid_coefficients( order=half_width )
    diagMat = np.zeros(shape=(N, N))
    for i, coef in enumerate(coefs):
        diagMat += np.diag([coef]*(N-abs(half_width-i)), half_width-i)
    smoothed_spectrum = np.dot(diagMat, spectrum)        
    return smoothed_spectrum

def Savitzky(spectrum: np.array, half_width: int=2, order: int=2):
    '''
    Savitzky-Golay smoothing method, return the smoothed spectrum
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
        -y(0)(fitted) = A(mth line) * M * y = (M * y)[0] 
        * M = (A.t*A).I *A.t, b = M * y
        
        -B = A * (A.t*A).I *A.t
        -y(fitted) = B * y = B * [y(m), y(m-1),...,y(-m)].T
        -y(0)(fitted) = B(mth line) * [y(m), y(m-1),...,y(-m)].T 
        * Above two solutions are coherent. 
        
    :param spectrum: analysing spectrum
    :param half_width: halfwidth of trandform window
    :param order: order of fitted polynomial
    :return smoothed: smoothed spectrum
    '''
    mat_order, mat_width = np.meshgrid(np.arange(order+1), np.arange(-half_width, half_width+1))
    A = mat_width ** mat_order
    M = np.matmul( np.linalg.inv(np.matmul(A.T, A)), A.T )
    smoothed_spectrum = spectrum.copy()
    for i in range(half_width, spectrum.shape[0]-half_width):
        smoothed_spectrum[i] = np.matmul( M, spectrum[i-half_width: i+half_width+1] )[0]
    return smoothed_spectrum

def Fourier(spectrum:np.array, type='low', fraction=0.2, sigma_erg=2, sigma_freq=0.01):
    '''
    Fourier smoothing using low-pass frequence filter
    
    :param spectrum: analysing spectrum
    :param type: type of low-pass function：'low', 'gauss', 'mixed'
    :param fraction: used for 'delta' and 'mixed', fraction of low frenquence perserved, between 0~0.5
    :param sigma_freq: used for 'gauss' and 'mixed', sigma of gauss distribution in frequence domain, between 0~0.1
    :param sigma_erg: used for 'gauss' and 'mixed', sigma of gauss distribution in energy domain, sigma_erg = len(spectrum) / simga_freq / 2pi 
    :return smoothed: smoothed spectrum
    '''
    
    if type == 'low':
        coefficients = np.ones(shape=spectrum.shape)
        coefficients[int(fraction*spectrum.shape[0]): int((1-fraction)*spectrum.shape[0])] = 0
    elif type == 'gauss':
        sigma = sigma_erg
        coefficients = np.arange(spectrum.shape[0]) 
        coefficients = coefficients**2 / spectrum.shape[0]**2 * -2 * np.pi**2 * sigma**2
        coefficients = np.exp(coefficients) + np.exp(coefficients)[::-1]
    elif type == 'mixed':
        sigma = sigma_erg
        coefficients = np.arange(spectrum.shape[0]) 
        coefficients = coefficients**2 / spectrum.shape[0]**2 * -2 * np.pi**2 * sigma**2
        coefficients = np.exp(coefficients) + np.exp(coefficients)[::-1]
        coefficients = coefficients / coefficients[int(fraction*spectrum.shape[0])]
        coefficients[0: int(fraction*spectrum.shape[0])] = 1
        coefficients[int((1-fraction)*spectrum.shape[0]): ] = 1
    else:
        coefficients = np.ones(shape=spectrum.shape)
    
    transformed_spectrum = np.fft.fft(spectrum)
    smoothed_transformed_spectrum = transformed_spectrum * coefficients
    smoothed_spectrum = np.real(np.fft.ifft( smoothed_transformed_spectrum ))
    return smoothed_spectrum

    
def Wavelet(spectrum, wavelet='sym8', level=8, mode='symmetric'):
    '''
    Wavelet smoothing using low-pass frequence filter
    
    :param spectrum: analysing spectrum
    :param type: type of wavelet function：
    :param order: order of wavelet tranform
    :param mode: padding mode
    :return smoothed: smoothed spectrum
    '''
    coefficients = pywt.wavedec(spectrum, wavelet=wavelet, mode=mode, level=level)
    sigma = np.median(np.abs(coefficients)[-1]) / 0.6745 # Donoho's rule
    threshold = sigma * (2*np.log(spectrum.shape[0])/np.log(10))**0.5 # VisuShrink
    # threshold = sigma * (2*np.log(spectrum.shape[0]))**0.5 # SureShrink
     
    for i in range(level-2, level):
        coefficients[i][abs(coefficients[i]) < threshold] = 0 # set high-order greater coefficients to oppisite shift [soft threshold]
        coefficients[i][coefficients[i] > +threshold] -= threshold 
        coefficients[i][coefficients[i] < -threshold] += threshold
        # coefficients[i][abs(coefficients[i]) > +threshold] = 0 # set high-order greater coefficients to zero[hard threshold]
        # coefficients[i][:] = 0 # set high-order coefficients to zero[hard threshold]
    smoothed_spectrum = pywt.waverec(coefficients,wavelet=wavelet,mode=mode)
    return smoothed_spectrum

def Translation_Invariance_Wavelet(spectrum, wavelet='sym8', level=8, mode='symmetric', step=1, times=None):
    '''
    Translation Invariance Wavelet smoothing using low-pass frequence filter 
    
    :param spectrum: analysing spectrum
    :param type: type of wavelet function：
    :param order: order of wavelet tranform
    :param mode: padding mode
    :param step: step of translation
    :param times: times of translation
    :return smoothed: smoothed spectrum
    '''
    if times is None:
        times = spectrum.shape[0]//step - 1

    spectrum = spectrum.copy()
    smoothed_spectrum_sum = np.zeros(shape=spectrum.shape[0])
    for i in range(1, times+1):
        spectrum2 = spectrum.copy()
        spectrum[:step], spectrum[step:] = spectrum2[-step:], spectrum2[:-step] # Translation
        smoothed_spectrum = Wavelet(spectrum, wavelet=wavelet, mode=mode, level=level) # Wavelet Reconstruction
        smoothed_spectrum2 = smoothed_spectrum.copy()
        step2 = (step*i) % spectrum.shape[0]
        smoothed_spectrum[-step2:], smoothed_spectrum[:-step2] = smoothed_spectrum2[:step2], smoothed_spectrum2[step2:] # Reverse Translation
        smoothed_spectrum_sum += smoothed_spectrum 
    smoothed_spectrum_sum /= (times) # Average the smoothed spectrum
    return smoothed_spectrum_sum

def SNIP(spectrum: np.array, half_width: int = 10):
    
    '''
    SNIP baseline stripping method

    :param spectrum: the analysing spectrum
    :param half_width：times of iteration also the halfwidth of transform window 

    :return baseline: the baseline of spectrum
    '''
    loglog_spectrum = np.log(np.log((spectrum+1)**0.5+1)+1)
    padded_spectrum = np.pad(loglog_spectrum, pad_width=half_width, mode='reflect')
    for p in range(half_width, 0, -1):
        padded_spectrum_compare = np.pad(padded_spectrum, pad_width=(p, 0), mode='edge')[: -p] / 2 + \
                                  np.pad(padded_spectrum, pad_width=(0, p), mode='edge')[ p: ] / 2 
        padded_spectrum = np.minimum(padded_spectrum, padded_spectrum_compare)
        
    baseline = ( np.exp( np.exp(padded_spectrum)-1 )-1 )**2 - 1
    baseline = baseline[half_width: -half_width]
    return baseline

def SNIP2(spectrum: np.array, info_calibration: tuple):

    '''
    SNIP baseline stripping method

    :param spectrum: the analysing spectrum
    :param half_width：times of iteration also the halfwidth of transform window 

    :return baseline: the baseline of spectrum
    '''

    loglog_spectrum = np.log(np.log((spectrum+1)**0.5+1)+1)
    padded_spectrum = np.pad(loglog_spectrum, pad_width=1, mode='reflect')
    half_width = int((info_calibration[0]*spectrum[-1]+info_calibration[1]) * 2.5)
    for p in range(half_width, 0, -1):
        left_index = int((p/2.5-info_calibration[1])/info_calibration[0])
        print(left_index)
        padded_spectrum_compare = np.pad(padded_spectrum, pad_width=(1, 0), mode='edge')[: -1] / 2 + \
                                  np.pad(padded_spectrum, pad_width=(0, 1), mode='edge')[ 1: ] / 2 
    padded_spectrum = np.minimum(padded_spectrum, padded_spectrum_compare)
    baseline = ( np.exp( np.exp(padded_spectrum)-1 )-1 )**2 - 1
    baseline = baseline[1: -1]
    return baseline


def evaluate(spectrum: np.array, smoothed_spectrum: np.array):
    power_signal  = np.mean(spectrum**2)
    power_noise = np.mean((smoothed_spectrum-spectrum)**2)
    SDR = 10*np.log10(power_signal/power_noise)
    RMSE = np.sqrt(np.mean((smoothed_spectrum-spectrum)**2))
    return SDR, RMSE


if __name__ == '__main__' :
    
    import matplotlib.pyplot as plt
    import simu_spectrum
    
    spectrum = simu_spectrum.simu_spectrum()
    plt.plot(spectrum)
    
    # baseline test
    # baseline = SNIP(spectrum, half_width=10)
    plt.plot(baseline)

    # for _ in range(1):
    #     spectrum1 = Wavelet(spectrum)
    #     spectrum2 = Translation_Invariance_Wavelet(spectrum)
    
    # print(evaluate(spectrum, spectrum1))
    # print(evaluate(spectrum, spectrum2))   
    # print(evaluate(spectrum, Centroid(spectrum))) 
    
    # plt.plot(spectrum1)
    # plt.plot(spectrum2)
    
    plt.show()