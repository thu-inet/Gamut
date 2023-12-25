def smooth(list_eff,width=11,order=4):
    '''
    Return smoothed spectrum with Savitzky-Golay filter

    :param list_eff: --mandatory, unsmoothed spectrum
    :param width: --Optional, window width, default 11, larger width leads to boardening
    :param order: --Optinal, Polynome order,default 4, higher order leads to stronger effect
    '''
    list_eff_smoothed=sgn.savgol_filter(list_eff,width,order)
    list_eff=list_eff_smoothed
    return list_eff
