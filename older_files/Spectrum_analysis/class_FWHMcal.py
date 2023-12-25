class FWHMcal():
    def __init__(self,form=-1,param=[]):
        self.form=form
        self.FWHM=self.FWHMfunction(param)
            
    def FWHMfunction(self,param):
        if self.form==0:
            return lambda erg: 1000*param[0]+param[1]*pow(erg*1000+param[1]*erg**2,0.5)
        elif self.form==-1:
            return lambda erg: 3
    # def FWHM(self,erg):
    #     a=-0.00139;b=0.00517;c=-0.376;erg=erg/1000
    #     s=max(erg+c*erg**2,0)       
    #     return (a+b*pow(s,0.5))*1000