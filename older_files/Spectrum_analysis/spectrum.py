# module name: GVpy

import class_ergcal
import class_FWHMcal

import func_import
import func_smooth
import func_search
import func_area
import func_others
import nuclib


class Spectrum():

    def __init__(self,list_erg=[],list_eff=[],list_eff_unsmoothed=[],ergcal=class_ergcal.Ergcal(0,0,0),fwhmcal=class_FWHMcal.FWHMcal(0,[-0.00139,0.00517,-0.376])):
        self.list_erg=list_erg
        self.list_eff=list_eff
        self.list_eff_unsmoothed=list_eff_unsmoothed
        self.ergcal=ergcal
        # self.fwhmcal=class_FWHMcal.FWHMcal(0,[0,0,0])
        self.fwhmcal=fwhmcal

    def import_GV(self,filepath): 
        list_erg,list_eff_unsmoothed,ergcal=func_import.import_GV(filepath,self.ergcal)
        self.list_erg=list_erg
        self.list_eff_unsmoothed=list_eff_unsmoothed
        self.ergcal=ergcal

    def import_MCNP(self,filepath):
        list_erg,list_eff_unsmoothed=func_import.import_MCNP(filepath,self.ergcal)
        self.list_erg=list_erg
        self.list_eff_unsmoothed=list_eff_unsmoothed
    
    def smooth_centroid(self,winlen=3):
        list_eff=func_smooth.smooth_centroid(self.list_eff_unsmoothed,winlen)
        self.list_eff=list_eff

    def smooth_fourier(self,type='low',factor=0.3):
        list_eff=func_smooth.smooth_fourier(self.list_eff_unsmoothed,type,factor)
        self.list_eff=list_eff

    def smooth_wavelet(self,wavelet='sym8',level=6,mode='symmetric'):
        list_eff=func_smooth.smooth_wavelet(self.list_eff_unsmoothed,wavelet,level,mode)
        self.list_eff=list_eff
    
    def SNIP(self):
        list_background=func_smooth.SNIP(self.list_eff_unsmoothed,ergcal=self.ergcal,fwhmcal=self.fwhmcal)
        self.list_background=list_background

    def add_peak(self,peak_left_boundary,peak_right_boundary):
        if not hasattr(self,'list_peaks'):
            self.list_peaks=[]
        self.list_peaks=func_search.add_peak(self.list_peaks,peak_left_boundary,peak_right_boundary)

    def search_simple(self,halfwidth=6):
        self.list_peaks=func_search.search_simple(self,halfwidth)

    def draw_spectrum(self,emin=0,emax=2400,ymin=1,ymax=1E5,NYlibpath=None):
        func_others.draw_spectrum(self,emin,emax,ymin,ymax,NYlibpath=NYlibpath)
    
    def strip(self,background,ratio):
        self.list_eff=func_others.strip(self,background,ratio)
    
    def areaTPA(self,erg,avg_number=4,width_base=3,base_method=0):
        net_counts_average, uncertainty, total_counts, peak_left_boundary, peak_right_boundary,base_start_left,base_start_right\
            =func_area.areaTPA(self,erg,avg_number=avg_number,width_base=width_base,base_method=base_method)
        return net_counts_average, uncertainty, total_counts, peak_left_boundary, peak_right_boundary,base_start_left,base_start_right
        
    def report(self,GVlibpath='D:\\User\\Reports\\Suspect.txt',NYlibpath=None):
        func_others.report(self,GVlibpath=GVlibpath,NYlibpath=NYlibpath)

    def export_GV(self,filepath,sample_discription='No sample description was entered.',live_time=10000,total_time=10000):
        func_others.export_GV(self,filepath,sample_discription=sample_discription,live_time=live_time,total_time=total_time)

    def match(self,GVlibpath='D:\\User\\Reports\\Suspect.txt',NYlibpath='E:\\Spectrum_Analysis\\Spectrum_data\\NYlib.xml.out'):
        gvlib=nuclib.GVlib(GVlibpath)
        self.list_gamma_match=gvlib.peak_match(list_peaks=self.list_peaks,ergcal=self.ergcal,NYlibpath=NYlibpath)

    def singlematch(self,energy,GVlibpath='D:\\User\\Reports\\Suspect.txt',NYlibpath='E:\\Spectrum_Analysis\\Spectrum_data\\NYlib.xml.out',tolerance=0.03,print_FOM=True):
        gvlib=nuclib.GVlib(GVlibpath)
        list_peaks=sorted(self.list_peaks,key=lambda peak: abs(self.ergcal.chn2erg(peak.centroid)-energy))[:1]
        gvlib.peak_match(list_peaks=list_peaks,ergcal=self.ergcal,NYlibpath=NYlibpath,print_FOM=print_FOM,tolerance=tolerance)

if __name__ == '__main__':
    Cs137=Spectrum()
    Cs137.import_GV('E:\\Spectrum_analysis\\Spectrum_data\\GM20燃耗球\\2020-07-21-1800s-023-目标球（反）.chn')
    Cs137.ergcal.recalibrate(list_chn=[2247,5058],list_erg=[661.66,1489.76],L=8192)
    Cs137.smooth_fourier(factor=0.3)
    Cs137.search_simple()
    # Cs137.singlematch(1384,NYlibpath='E:\\NUIT_data\\input\\inp_test.xml.out')
    Cs137.match(NYlibpath='E:\\NUIT_data\\input\\inp_HJ.xml.out')
    # Cs137.draw_spectrum(NYlibpath='E:\\NUIT_data\\input\\inp_test.xml.out')
    Cs137.report()
    x=1