import matplotlib.pyplot as plt
import numpy as np

def SNIP(list_erg: np.array, list_eff: np.array, M: int=10):
    
    '''
    Background substration with SNIP method


    :param list_erg: --mandatory, energy list 
    :param list_eff: --mandatory, counts per channel list of same length with list_erg

    :return: spectrum of background in length of list_eff
    '''
    LLSs = np.log( np.log( list_eff + np.ones(list_eff.shape[0])) + np.ones(list_eff.shape[0]) )
    LLSs = np.concatenate( (LLSs[:M], LLSs, LLSs[list_eff.shape[0]-M: ]) )
    L = LLSs.shape[0]
    for j in range(M):
        LLSs_new = (LLSs[0: L-2*j] + LLSs[2*j: L]) / 2 
        updates = np.where( ( LLSs_new - LLSs[j: L-j]) <0 )[0]
        LLSs[updates + np.full(updates.shape, j)] = LLSs_new[updates]
    smoothed  = np.exp( np.exp(LLSs) - 1 ) -1
    smoothed = smoothed[M: L-M]
    
    return smoothed



if __name__=='__main__':
    
    
    import matplotlib.pyplot as plt
    f = lambda i, j: np.exp(-i**2/20)*j
    list_erg = np.linspace(0, 1000, 500)
    list_eff = np.random.normal(loc=1000, scale=40, size=(500,)) \
        + np.array([f(i-300, 1000) for i in range(500)]) \
        + np.array([f(i-100, 300) for i in range(500)])  \
        + np.array([f(i-115, 500) for i in range(500)])  \
        + np.array([f(i-415, 200) for i in range(500)])
    plt.plot(list_erg, list_eff)
    # SNIP(np.arange(0,10,1), np.array([1,4,2,6,5,7,7,4,3,2]), M=3)
    list_eff = SNIP(list_erg, list_eff, M=20)
    plt.plot(list_erg, list_eff)
    plt.show()

