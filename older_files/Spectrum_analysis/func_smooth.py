import numpy as np

def recur(i: int):
    '''
    Auxiliary funciton, calculate centroid smoothing coefficients by recursion.
    
    :param i: index, exit the recurtion at i=1
    :return mat: matrix 
    '''
    if i == 0:
        return np.array([1])
    else:
        mat = recur(i-1)
        mat1 = np.concatenate( (mat, np.zeros(2)) )
        mat2 = np.concatenate( (np.zeros(1), mat, np.zeros(1)) )
        mat3 = np.concatenate( (np.zeros(2), mat) )
        return mat1/4 + mat2/2 + mat3/4
    
    
def centroid(list_eff: np.array, L: int = 3, T: int =1):
    '''
    centroid smoothing method, return smoothed spectrum
    
    :param list_eff: Mandatory, the analysing spectrum
    :param L: Optional, the halfwidth of transform window 
    :param T: Optional, time to repeat the smoothing
    
    :return smoothed: smoothed spectrum
    '''
    trans = recur(L) 
    for _ in range( T ):   
        transed = np.array( [ np.dot( trans, list_eff[i-L: i+L+1] ) for i in range( L, len(list_eff)-L ) ] )
        list_eff = np.concatenate( (transed[0: L], transed, transed[-L: ]) ) # extrapolate the borders
    return list_eff

def savol(list_eff: np.array, L: int =4, R: int=2):
    '''
    Savitzy-Golay smoothing method, return the smoothed spectrum
    Equation: 
        Ab+E=x
        -b=[b0,b1,b2,...,bn], n=order of polymomial, b1,...,bn=polynomial coefficients
        -A=[[1   m     m^2    ,...,       m^n], 
            [1  m-1  (m-1)^2  ,...,   (m-1)^n]
            [1   i     i^2    ,...,       i^n]
            [1  -m   (-m)^2   ,...,    (-m)^n]], m=half width of transformation window
            * This format
        -x=[y(m),y(m-1),...,y(0),...,y(-m)]
    Solution: 
        -b=(At*A)^(-1)*At*x
        -x=Ab=A*(At*A)^(-1)*At*x

        -B=A*(At*A)^(-1)*At, x(new)=B*x=B*[y(m),y(m-1),...,y(-m)]t
        -xi(new)=B(ith row)*x

    '''
    A = np.array( [ [ i**j for j in range(R+1) ] for i in range(-L, L+1) ] )
    b = np.matmul( np.linalg.inv(np.matmul(A.T, A)), A.T )
    smoothed = list_eff.copy()
    for i in range(L, list_eff.shape[0]-L):
        smoothed[i] = np.matmul( b, list_eff[i-L:i+L+1] )[0]
    return smoothed

def fourier(list_eff_unsmoothed:np.array, type='low', f=0.3, sigma=3):

    transformed = np.fft.fft(list_eff_unsmoothed)
    smoothed_transformed = transformed.copy()
    if type == 'low':
        l = int( len(list_eff_unsmoothed)/2 * f )
        for i in range(l, len(smoothed_transformed)-l):
            smoothed_transformed[i] = 0
    elif type == 'gauss':
        k = 2 * np.pi**2 / len(smoothed_transformed)**2 * sigma**2
        for i in range(len(smoothed_transformed)):
            j = int(len(smoothed_transformed)/2) - int(abs(len(smoothed_transformed)/2-i))
            smoothed_transformed[i] *= np.exp( -k * j**2 )
    elif type == 'mixed':
        l = int( len(list_eff_unsmoothed)/2 * f )
        k = 2 * np.pi**2 / len(smoothed_transformed)**2 * sigma**2
        for i in range(len(smoothed_transformed)):
            if ( i>l ) and ( i<len(smoothed_transformed)-l ):
                j = int(len(smoothed_transformed)/2) - int(abs(len(smoothed_transformed)/2-i))
                smoothed_transformed[i] *= np.exp( -k * j**2 )
    smoothed = np.fft.ifft( smoothed_transformed )
    smoothed = np.real(smoothed)
    return smoothed

# def CVFTM(list_eff_unsmoothed):
#     list_odd = np.array([list_eff[i] for i in range(1, len(list_eff_unsmoothed, 2))])
#     list_even = np.array([list_eff[i] for i in range(0, len(list_eff_unsmoothed, 2))])
#     interpolated_odd = list_odd[1:] + list_odd[:-1]
#     interpolated_even = list_odd[1:] + list_odd[:-1]
#     transformed_odd = np.fft.fft(list_odd)
#     transformed_even = np.fft.fft(list_even)
#     for f in np.linspace(0, 0.25, 0.01):
        
        # smoothed_transformed = transformed.copy()
        # if type=='low':
        #     l = int( len(list_eff_unsmoothed)/2*(1-f) )
        #     for i in range(l, len(smoothed_transformed)-l):
        #         smoothed_transformed[i]=0
        # elif type=='gauss':
        #     f = f / len(smoothed_transformed) / 2
        #     for i in range(len(smoothed_transformed)):
        #         smoothed_transformed[i] /= np.exp( f * i**2 )
        # smoothed = np.fft.ifft( smoothed_transformed )
        # smoothed = np.real(smoothed)
        # smoothed_even = fourier(list_even, f)[:-1]
        # smoothed_odd = fourier(list_odd, f)[1:]
        # error = (interpolated_odd-smoothed_odd)**2 * (interpolated_even-smoothed_even)**2
    
    
    

def smooth_wavelet(list_eff_unsmoothed,wavelet='sym8',level=6,mode='symmetric'):

    coeffs=wavedec(list_eff_unsmoothed,wavelet=wavelet,mode=mode,level=level)
    sigma=np.median([abs(s) for s in coeffs[-1]])/0.6745*(np.log(len(list_eff_unsmoothed))*2)**0.5
    for i_cDl in range(1,len(coeffs)):
        for i_cD in range(len(coeffs[i_cDl])):
            detail=coeffs[i_cDl][i_cD]
            if detail > sigma:
                coeffs[i_cDl][i_cD]-=sigma
            elif detail < -sigma:
                coeffs[i_cDl][i_cD]+=sigma
            else:
                pass
    list_eff=waverec(coeffs,wavelet=wavelet,mode=mode)[:-1]
    return list_eff


def SNIP(list_eff: np.array, M: int = 10):
    
    '''
    SNIP baselien substraction method, return the baseline spectrum

    :param list_eff: the analysing spectrum
    :param M：Optional, times of iteration also the halfwidth of transform window 

    :return: spectrum baseline
    '''
    
    LLSs = np.log( np.log( list_eff + np.ones(list_eff.shape[0])) + np.ones(list_eff.shape[0]) )
    LLSs = np.concatenate( (LLSs[:M], LLSs, LLSs[list_eff.shape[0]-M: ]) )
    L = LLSs.shape[0]
    
    for j in range(M):
        LLSs_new = ( LLSs[0: L-2*j] + LLSs[2*j: L] ) / 2 
        updates = np.where( LLSs_new < LLSs[j: L-j] )[0]
        LLSs[updates + np.full(updates.shape, j)] = LLSs_new[updates]
    baseline  = np.exp( np.exp(LLSs) - 1 ) -1
    baseline = baseline[M: L-M]
    
    return baseline

if __name__ == '__main__' :
    
    import simu_spectrum2
    import matplotlib.pyplot as plt
    
    f = lambda i, j: np.exp(-i**2/20)*j
    list_erg = np.linspace(0, 1000, 500)
    list_eff = np.random.normal(loc=1000, scale=40, size=(500,)) \
        + np.array([f(i-300, 1000) for i in range(500)]) \
        + np.array([f(i-100, 300) for i in range(500)])  \
        + np.array([f(i-115, 500) for i in range(500)])  \
        + np.array([f(i-415, 200) for i in range(500)])
    plt.plot(list_erg, list_eff)
    
    list_eff = centroid(list_eff, L=5, T=3)
    
    plt.plot(list_erg, list_eff)
    plt.show()