from numpy import polyfit

class Ergcal():
    def __init__(self,c0,c1,c2=0):
        self.c0=c0
        self.c1=c1
        self.c2=c2
        if c2==0:
            self.type='linear' 
        else:
            self.type=='quadratic'

    def chn2erg(self,chn):
        return self.c2*chn**2+self.c1*chn+self.c0
    
    def erg2chn(self,erg):
        if self.type=='linear':
            return int((erg-self.c0)/self.c1)
        else:
            return int((self.c1-(self.c1**2-4*self.c2*(self.c0-erg))**0.5)/2/self.c2)

    def recalibrate(self,L,list_chn,list_erg,type='linear'):

        self.type=type
        if len(list_erg)!=len(list_chn):
            raise ValueError('ValueError --Lengths don\'t match!')
        if type=='linear':
            param=polyfit(list_chn,list_erg,deg=1)
            self.c2=0
            self.c1=param[0]
            self.c0=param[1]
        else:
            param=polyfit(list_chn,list_erg,deg=2)
            self.c2=param[0]
            self.c1=param[1]
            self.c0=param[2]
        list_erg_new=[self.chn2erg(chn) for chn in range(L)]
        return list_erg_new


if __name__=='__main__':

    # Initiate ergcal instance
    ergcal=Ergcal(c0=-5.391662E-001,c1=2.947627E-001,c2=-4.814516E-008)
    list_erg=[ergcal.chn2erg(chn) for chn in [100,1000,4000,8000]]
    print('Quadratic:',list_erg)
    '''
    Quadratic: [28.9366223484, 294.17538864, 1177.74131124, 2354.48114356]
    '''

    # Recalibrate with chn-erg lists Order 2 
    list_chnx=[2247.21,2600.99,7425.82]
    list_ergx=[661.61,765.81,2185.66]
    ergcal.recalibrate(list_chnx,list_ergx,type='linear')
    list_erg=[ergcal.chn2erg(chn) for chn in [100,1000,4000,8000]]
    print('Linear coefficients:',ergcal.c0,ergcal.c1)
    print('Recalibrate, linear1:',list_erg)
    '''
    Linear coefficients: 0.32086519719396356 0.294289678883399
    Recalibrate, linear1: [29.749833085533865, 294.61054408059294, 1177.47958073079, 2354.638296264386] 
    '''

    # Recalibrate with chn-erg lists Order 3
    ergcal.recalibrate(list_chnx,list_ergx,type='quadratic')
    list_erg=[ergcal.chn2erg(chn) for chn in [100,1000,4000,8000]]
    print('Quadratic coefficients:',ergcal.c0,ergcal.c1,ergcal.c2)
    print('Recalibrate, quadratic:',list_erg)
    '''
    Quadratic coefficients: -0.5543851793937203 0.29477067736930873 -4.895663826162909e-08
    Recalibrate, quadratic: [28.92219299115454, 294.16733555165337, 1177.7450180856551, 2354.477808926332]
    '''
    


