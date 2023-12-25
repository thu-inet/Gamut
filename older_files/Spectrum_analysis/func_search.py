import class_objects
import numpy as np
from scipy.signal import find_peaks, argrelextrema

def add_peak(listpeak,peakleft_boundary,peak_right_boundary):
    centroid=int((peakleft_boundary+peak_right_boundary)/2)
    peak=class_objects.Peak(centroid,peakleft_boundary,peak_right_boundary,gamma='')
    listpeak.append(peak)
    return listpeak

def search_simple(list_eff: np.array, halfwidth: int = 6, list_peaks: list=[]):

    IF = np.sign(-list_eff[0:-1] + list_eff[1:])  # calculate IF-increment factor
    IF = np.concatenate((IF[:1], IF))
    IFX = np.sign(0.6*IF[0:-2] + IF[1:-1] + 0.6*IF[2:]) # correct outliers
    IFX = np.concatenate((IFX[:1], IFX, IFX[-1:])) 
    IFY = [((sum(IFX[i-halfwidth:i]) == halfwidth) and (sum(IFX[i:i+halfwidth]) == -halfwidth)).astype(int) for i in range(halfwidth, IFX.shape[0]-halfwidth)] 
    IFY = np.concatenate((IFY[:halfwidth],IFY, IFY[-halfwidth:])) # calculate bool: consecutive increassing points and decreasing points
    IFZ = np.where(IFY==1)[0] - 1 # -1 to correct displacement
    for centroid in IFZ:
        try:
            left_len = min(centroid, 3*halfwidth) - np.where(IFX[max(centroid-3*halfwidth,0): centroid] == -1)[0][-1]
        except:
            left_len = 3 * halfwidth
        try:
            right_len = np.where(IFX[centroid: centroid+3*halfwidth] == 1)[0][1]
        except:
            right_len = 3 * halfwidth
        peak = (max(centroid-left_len, 0),
             centroid, 
             min(centroid+right_len, list_eff.shape[0]))
        list_peaks.append(peak)
    return list_peaks
    

def gauss(list_eff: np.array, tol: float = 7, order: int=2, list_peaks: list = []):

    values = list_eff[2:-order] * list_eff[order:-2]  / list_eff[0:-order-2] / list_eff[2+order :] # calcultate Gauss multiplication value
    values = (values - np.ones(values.shape)) * list_eff[2:-order]**0.5 # calculate the prominence
    values = np.concatenate((values[:2], values))
    indexes = np.where(values > tol)[0]  # filter out all values that exceed tolerance, which means peaks
    splits = np.where(indexes[1:]-indexes[:-1]>2)[0] # filter out all indexes that are not consecutive, which means seperate peaks
    splits = np.concatenate((splits, np.array([-1]))) # indexes=[...peakleft_boundary, ..., peak_right_boundary[index=split], peakleft_boundary[index=split+1], ..., ]
     
    left = indexes[0]
    for i in range(len(splits)) : # first peak region, starts with indexes[0]
        right = indexes[splits[i]] # split means splitting indexes of indexes
        centroid = int((left + right) / 2) # centroid of the peak region by middle point of left, _right
        if values[centroid] > tol/2  : 
            sigma = np.log( values[centroid] / list_eff[centroid]**0.5+1 )**(-0.5) # sigam of the centroid
            peak = (
                max(int(centroid-3*sigma), 0), 
                centroid, min(int(centroid+3*sigma), 
                list_eff.shape[0]))
            list_peaks.append(peak)
        left = indexes[splits[i]+1]

    return list_peaks

def derive(list_eff: np.array, mode: int = 322, list_peaks: list = []):

    '''
    Numerical differential method. Returns peak list in the spectrum

    :param list_eff: Compulsary, the analysing spectrum
    :param mode: Optional, integer with 3/4 digits: [Order of numerical differential|Order of fitting polunomial|halfwidth of convolutional window|Repeated times]
    :param list_peaks: Optional, necessary when SZAC is called for adding peaks

    :output list_peaks: peak list, elements: tuple (left, centroid, right) 
    '''
    list_eff = np.array(list_eff)
    if mode == 122:
        diff = -1/5*list_eff[0:-4] - 1/10*list_eff[1:-3] + 1/10*list_eff[3:-1] + 1/5*list_eff[4:] # diffS1N2M2
        diff = np.concatenate((diff[:2], diff))
        diff = np.concatenate((diff, diff[-2:]))
    elif mode == 322:
        diff = -0.5*list_eff[0:-4] - list_eff[1:-3] + list_eff[3:-1] + 0.5*list_eff[4:] # diffS3N2M2
        diff = np.concatenate((diff[:2], diff))
        diff = np.concatenate((diff, diff[-2:]))
    elif mode == 1325:
        def S1N3M2(list_eff): # diffS1N3M2 连续使用三次一阶三次五点光滑寻峰
            diff = 1/12*list_eff[0:-4] - 8/12*list_eff[1:-3] + 8/12*list_eff[3:-1] - 1/12*list_eff[4:] 
            diff = np.concatenate((diff[:2], diff))
            diff = np.concatenate((diff, diff[-2:]))
            return diff
        diff1 = S1N3M2(list_eff)
        diff2 = S1N3M2(diff1)
        diff3 = S1N3M2(diff2)
        diff = diff3
    else:
        diff = -0.5*list_eff[0:-4] - list_eff[1:-3] - list_eff[3:-1] + 0.5*list_eff[4:] # diffS3N2M2
        diff = np.concatenate((diff[:2], diff))
        diff = np.concatenate((diff, diff[-2:]))

    plt.plot(diff,'-x', color='black')
    plt.show()

    # indexes = np.where(diff > 10)[0]
    # splits = np.where((indexes[1:]-indexes[:-1]) >2)[0]
    # splits = np.concatenate((splits, np.array([-1])))
    # left = indexes[0]
    # for i in range(splits.shape[0]):
    #     right = indexes[splits[i]]
    #     centroid = int(0.5*(left+right))
    #     peak = (centroid-H*2, centroid, centroid+H*2)
    #     list_peaks.append(peak)
    #     left = indexes[splits[i]+1]

    return list_peaks
    

def SZAC(list_eff: np.array, H: int=4, m: int=5, f: float=5, list_peaks: list=[], func: str = 'Gauss'):

    '''
    Sysmetric Zero Area Convolution Method. Returns peak list in the spectrum

    :param list_eff: Compulsary, the analysing spectrum
    :param H: Optional, FWHM of peaks, applied in built-in function constructions
    :param m: Optional, haldwidth of convolution window, big m(5) to suppress the baseline, small m(1) to distinguish overlapping peaks
    :param list_peaks: Optional, necessary when SZAC is called for adding peaks
    :param func: Optional, allow user defined peakshape-like functions, default is Gauss, Cauchy and cos2 supported

    :output list_peaks: peak list, elements: tuple (left, centroid, right)      
    '''

    if func == 'Gauss':
        func = lambda x: np.exp(-(x/H)**2) * 1/16 
    elif func == 'cos2':
        func = lambda x: np.cos(x/H * np.pi/2)**2
    elif func == 'Cauchy':
        func = lambda x: 1/(1+(x/H)**2)
    elif func == 'sech':
        func = lambda x: 1/np.cosh(2.634*x/H)
    else:
        pass

    d = sum([func(j) for j in range(-m, m+1)]) / (2*m +1)
    mat1 = np.array([ [ func((i-j)) - d     if abs(j-i) <=m else 0 for j in range(list_eff.shape[0])] for i in  range(list_eff.shape[0])])
    mat2 = np.array([ [(func((i-j)) - d)**2 if abs(j-i) <=m else 0 for j in range(list_eff.shape[0])] for i in  range(list_eff.shape[0])])
    vals1 = np.dot(mat1, list_eff)
    vals2 = np.dot(mat2, list_eff)
    vals1[:m], vals1[-m:] = vals1[m:2*m], vals1[-2*m:-m]
    vals2[:m], vals2[-m:] = vals2[m:2*m], vals2[-2*m:-m]
    SS = vals1 / vals2**0.5

    indexes = np.where(SS > f)[0]
    splits = np.where((indexes[1:]-indexes[:-1]) >2)[0]
    splits = np.concatenate((splits, np.array([-1])))
    left = indexes[0]
    for i in range(splits.shape[0]):
        right = indexes[splits[i]]
        centroid = int(0.5*(left+right))
        peak = (centroid-H*2, centroid, centroid+H*2)
        list_peaks.append(peak)
        left = indexes[splits[i]+1]

    return list_peaks

def doubleSZAC(list_eff: np.array):
    
    '''
    Apply SZAC two times, first with m=5 to suppress the baseline and second with m=1 to distinguish overlapping peaks 
    '''
    list_peaks = SZAC(list_eff, m=5, func='Gauss')
    list_peaks_new = []
    for peak in list_peaks:
        left, centroid, right = peak
        list_peaks_new = SZAC(list_eff[left: right], m=1, list_peaks= list_peaks_new)
    return list_peaks_new


if __name__ == '__main__':

    import matplotlib.pyplot as plt
    f = lambda i, j: np.exp(-i**2/20)*j
    list_eff = np.random.normal(loc=1000, scale=40, size=(500,)) \
        + np.array([f(i-300, 1000) for i in range(500)]) \
        + np.array([f(i-100, 300) for i in range(500)])  \
        + np.array([f(i-115, 500) for i in range(500)])  \
        + np.array([f(i-415, 200) for i in range(500)])
    # list_peaks = search_simple(list_eff)
    # list_peaks = gauss(list_eff)
    list_peaks = derive(list_eff)
    # list_peaks = SZAC(list_eff, func='Cauchy')
    plt.plot(list_eff,linestyle='-', color='black')
    for peak in list_peaks:
        left, centroid, right = peak
        plt.plot(np.arange(left, min(right+1, list_eff.shape[0])), list_eff[left: right+1], 'x')
    plt.savefig(r'C:\Users\Alber\Desktop\x.jpg')
    plt.show()
    