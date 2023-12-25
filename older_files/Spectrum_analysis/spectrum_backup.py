# module name: GVpy
from numpy import polyfit
import matplotlib.pyplot as plt
from scipy.fftpack import fft, ifft
import numpy as np
from pywt import wavedec, waverec
import time
import scipy.signal as sgn

import class_ergcal
import class_FWHMcal
import class_objects

class Spectrum():

    def __init__(self):
        self.list_erg=[]
        self.list_eff=[]
        self.list_eff_unsmoothed=[]
        self.ergcal=class_ergcal.Ergcal(0,0,0)
        self.fwhmcal=class_FWHMcal.FWHMcal(0,[-0.00139,0.00517,-0.376])

    def import_GV(self,filepath): 

        self.list_erg=[]
        self.list_eff=[]
        self.list_eff_unsmoothed=[]

        if not filepath[-3:].lower() in ['spe','chn']:
            raise ValueError('Error-This file cannot be interpreted as Spe or MCNP-ouput format.')
        else:
            pass

        with open(filepath,'r',encoding='utf-8') as fileopen:
            filelines=fileopen.readlines()

        indl=0
        while True:
            line=filelines[indl]
            if '$DATA' in line:
                channels=int(filelines[indl+1].split()[1])
                indl+=2
                break
            else:
                indl+=1

        for i in range(indl,indl+channels):
            line=filelines[i]
            self.list_eff_unsmoothed.append(int(line[:-1]))
            self.list_eff.append(int(line[:-1]))
        indl+=channels

        while True:
            line=filelines[indl]
            if '$ENER_FIT' in line:
                indl+=1
                break
            else:
                indl+=1
        
        line=filelines[indl]
        self.ergcal.c0=float(line.split()[0])
        self.ergcal.c1=float(line.split()[1])
        self.list_erg=[self.ergcal.chn2erg(chn) for chn in range(len(self.list_eff))]

        return 0

    def import_MCNP(self,filepath):

        self.list_erg=[]
        self.list_eff=[]
        self.list_eff_unsmoothed=[]

        with open(filepath,'r',encoding='utf-8') as fileopen:
            filelines=fileopen.readlines()

        while True:
            line=filelines[indl]
            if 'tally   8' in line:
                indl+=6
                break
            else:
                indl+1

        while True:
            line=filelines[indl]
            try:
                erg,eff,err=line.split()
                self.list_erg.append(float(erg))
                self.list_eff_unsmoothed.append(float(eff))
                self.list_eff.append(float(eff))
                indl+=1
            except:
                break
        self.list_eff_unsmoothed=self.list_eff
        return 0
    
    # def FWHM(self,erg):
    #     a=-0.00139;b=0.00517;c=-0.376;erg=erg/1000
    #     s=max(erg+c*erg**2,0)       
    #     return (a+b*pow(s,0.5))*1000

    def smooth_centroid(self,winlen=3):
        list_eff_unsmoothed=self.list_eff_unsmoothed
        if winlen==3:
            for i in range(1,len(self.list_eff)-1):
                self.list_eff[i]=(list_eff_unsmoothed[i-1]+2*list_eff_unsmoothed[i]+list_eff_unsmoothed[i+1])/4
        elif winlen==5:
            for i in range(2,len(self.list_eff)-2):
                self.list_eff[i]=(list_eff_unsmoothed[i-2]+4*list_eff_unsmoothed[i-1]+6*list_eff_unsmoothed[i]+4*list_eff_unsmoothed[i+1]+list_eff_unsmoothed[i+2])/16
        elif winlen==7:
            for i in range(3,len(self.list_eff)-3):
                self.list_eff[i]=(list_eff_unsmoothed[i-3]+6*list_eff_unsmoothed[i-2]+15*list_eff_unsmoothed[i-1]+20*list_eff_unsmoothed[i]+15*list_eff_unsmoothed[i+1]+6*list_eff_unsmoothed[i+2]+list_eff_unsmoothed[i+3])/64
        else:
            for i in range(1,len(self.list_eff)-1):
                list_eff_unsmoothed[i]=(list_eff_unsmoothed[i-1]+2*list_eff_unsmoothed[i]+list_eff_unsmoothed[i+1])/4
        return 0

    def smooth_fourier(self,type='low',factor=0.3):
        trsed=fft(self.list_eff_unsmoothed)
        trsed_smthed=trsed
        if type=='low':
            l=int(len(self.list_eff_unsmoothed)/2*(1-factor))
            for i in range(l,len(trsed_smthed)-l):
                trsed_smthed[i]=0
        elif type=='gauss':
            factor=factor/len(trsed_smthed)/2
            for i in range(len(trsed_smthed)):
                trsed_smthed[i]*=np.exp(-factor*i**2)
        smthed=ifft(trsed_smthed)
        self.list_eff=np.real(smthed)
        return 0

    def smooth_wavelet(self,wavelet='sym8',level=6,mode='symmetric'):
        coeffs=wavedec(self.list_eff_unsmoothed,wavelet=wavelet,mode=mode,level=level)
        sigma=np.median([abs(s) for s in coeffs[-1]])/0.6745*(np.log(len(self.list_eff_unsmoothed))*2)**0.5
        for i_cDl in range(1,len(coeffs)):
            for i_cD in range(len(coeffs[i_cDl])):
                detail=coeffs[i_cDl][i_cD]
                if detail > sigma:
                    coeffs[i_cDl][i_cD]-=sigma
                elif detail < -sigma:
                    coeffs[i_cDl][i_cD]+=sigma
                else:
                    pass
        self.list_eff=waverec(coeffs,wavelet=wavelet,mode=mode)[:-1]
        return 0
    

    def SNIP(self):
        list_v=[0]*len(self.list_erg)
        for i in range(len(list_v)):
            list_v[i]=np.log(np.log((self.list_eff_unsmoothed[i]+1)**0.5+1)+1)

        centroid=100
        width=100
        while True:
            if hasattr(self,'fwhmcal'):
                m=max(int(self.fwhmcal(self.ergcal.chn2erg(centroid))/self.ergcal.c1/2.35*3.5),2)
            else:
                m=10
            for p in range(1,m+1):
                if centroid==width:
                    for i in range(p,centroid+width+1):
                        list_v[i]=min(list_v[i],(list_v[i-p]+list_v[i+p])/2)
                elif len(list_v)-centroid<=width:
                    for i in range(centroid-width,len(list_v)-p):
                        list_v[i]=min(list_v[i],(list_v[i-p]+list_v[i+p])/2)
                else:
                    for i in range(centroid-width,centroid+width+1):
                        list_v[i]=min(list_v[i],(list_v[i-p]+list_v[i+p])/2)
            if len(list_v)-centroid<=width:
                break
            centroid+=2*width
        self.list_background=[]
        for i in range(len(list_v)):
            self.list_background.append((np.exp(np.exp(list_v[i])-1))**2-1)
        return 0


    def modify_peak(self,peak_left_boundary,peak_right_boundary):
        if not hasattr(self,'list_peaks'):
            self.list_peaks=[]
        centroid=int((peak_left_boundary+peak_right_boundary)/2)
        peak=class_objects.Peak(centroid,peak_left_boundary,peak_right_boundary,gamma='')
        self.list_peaks.append(peak)
        return 0

    def search_simple(self,halfwidth=4):
        list_IF=[1 if self.list_eff[i]<self.list_eff[i+1] else -1 for i in range(len(self.list_eff)-1)]
        list_IF.append(1)
        list_IFX=[(list_IF[i-1]+list_IF[i]+list_IF[i+1])/abs(list_IF[i-1]+list_IF[i]+list_IF[i+1]) for i in range(1,len(self.list_eff)-1)]
        list_IFX.insert(0,1)
        list_IFX.append(1)
        list_IF=list_IFX
        i=0
        if not hasattr(self,'list_peaks'):
            self.list_peaks=[]
        while True:
            try:
                if list_IF[i]==1:
                    j=1
                    while list_IF[i+j]==1:
                        j+=1
                    k=1
                    while list_IF[i+j+k]==-1:
                        k+=1
                    if  (j>=halfwidth) & (k>=halfwidth) & (self.list_eff[i+j]>self.list_eff[i]+0.5*max(self.list_eff[i],0)**0.5):
                        peak=class_objects.Peak(centroid=i+j,peak_left_boundary=i,peak_right_boundary=min(i+j+k,len(self.list_eff)),gamma='')
                        self.list_peaks.append(peak)
                    i+=j+k+1
                else:
                    i+=1
            except:
                return 0
    
    # def search_SZAC(self,halfwindowwidth=8):
    #     list_IF=[1]*len(self.list_eff)
    #     for i in range(halfwindowwidth,len(self.list_eff)-halfwindowwidth):
    #         IF=2*self.list_eff[i]
    #         for j in range(-halfwindowwidth,halfwindowwidth+1):
    #             IF-=self.list_eff[i+j]
    #         list_IF[i]=IF
    #     return list_IF


    # def gauss(list_eff,order=1,factor=1):
    #     list_PG=[list_eff[i]*list_eff[i+order-2]/list_eff[i-2]/list_eff[i+order]-1 for i in range(2,len(list_eff)-order)]
    #     list_PG.insert(0,0)
    #     list_PG.insert(1,0)
    #     for i in range(order):
    #         list_PG.append(0)
    #     list_peaks=[]
    #     for i in range(len(list_PG)):
    #         if list_PG[i]>=factor/max(list_eff[i],1)**0.5:
    #             r=0
    #             while 1:
    #                 if (list_PG[i+r]+factor/list_eff[i])*(list_PG[i+r+1]+factor/list_eff[i])<=0:
    #                     break
    #                 else:
    #                     r+=1
    #             l=0
    #             while 1:
    #                 if (list_PG[i-l-1]+factor/list_eff[i])*(list_PG[i-l]+factor/list_eff[i])<=0:
    #                     break
    #                 else:
    #                     l-=1
    #             peak=Peak(centroid=i,left=i-l-1,right=i+r,gamma='')
    #             list_peaks.append(peak)
    #     return list_peaks

    def draw_spectrum(self,emin=0,emax=2400,ymin=1,ymax=10000):

        plt.figure(figsize=(30,10))
        plt.plot(self.list_erg,self.list_eff_unsmoothed,'b.',label='Unsmoothed',markersize=0.5)
        plt.plot(self.list_erg,self.list_eff,'b-',label='Smoothed',linewidth=0.5)
        for peak in self.list_peaks:
            centroid=list(self.list_eff[peak.peak_left_boundary:peak.peak_right_boundary+1]).index(max(self.list_eff[peak.peak_left_boundary:peak.peak_right_boundary+1]))+peak.peak_left_boundary
            if self.list_eff[centroid]-self.list_eff[peak.peak_left_boundary]>20:
                plt.text(self.list_erg[centroid],self.list_eff_unsmoothed[centroid],'  <--%6.2fkeV(???)'%(self.list_erg[peak.centroid]),rotation=90,ha='center',va='bottom')
            plt.plot(self.list_erg[peak.peak_left_boundary:peak.peak_right_boundary+1],self.list_eff[peak.peak_left_boundary:peak.peak_right_boundary+1],'r-',linewidth=0.7)
            plt.fill_between(self.list_erg[peak.peak_left_boundary:peak.peak_right_boundary+1],self.list_eff[peak.peak_left_boundary:peak.peak_right_boundary+1],alpha=0.7)
        # plt.axvline(x=list_erg[peak_boundary_left],ymax=list_eff[peak_boundary_left]/ymax,linestyle='--')
        # plt.axvline(x=list_erg[peak_boundary_right],ymax=list_eff[peak_boundary_right]/ymax,linestyle='--')
        plt.xlabel('Energy, keV')
        plt.ylabel('Counts')
        plt.yscale('log')
        plt.title('Peak Region Plot')
        plt.xticks(np.arange(0,2500,100))
        plt.legend()
        plt.xlim([emin,emax])
        ymax=np.exp(int(np.log(max(self.list_eff))+1))
        plt.ylim([ymin,ymax])
        plt.savefig('C:\\Users\\Alber\\Desktop\\WJ.png')
        plt.show()
    
    def strip(self,background,ratio):
        list_eff_strpied=[]
        for total,back in zip(self.list_eff,background.list_eff):
            if total !=0:
                net=max(int(total-back*ratio+1),0)  
            else:
                net=0
            list_eff_strpied.append(net)
        self.list_eff_unsmoothed=list_eff_strpied
        self.list_eff=list_eff_strpied
        return 0
    
    
    def areaTPA(self,erg,avg_number=4,width_base=3,base_method=0):

        '''
        Avaraged total peak area method.
        Total peak area method: two sets of boundaries: one for base signal determination, one for peak area calculation
        1.Detemine the boundary_left and boundary_right by minimum or counts difference;
        2.avg_number, set by user, to determine how many times one should calculate the peak area by using a different baseline each time
        If avg_number equals 10, then the base will be determined using peak_boundary_left-1, -2,-3,...,-10 and peak_boundary_right+1,+2,...,+10.


        :param list_erg: --mandatory, energy list 
        :param list_eff: --mandatory, list of counts per channel OR list of effience per channel. length is preferred to be 8192
        :param peak: --mandatory, channel of the peak, slight error can be acccpted
        :param avg_number: --optional, times to repeat peak area calculation, higher avg_number enhances the result's validity but enlarges the uncertainty as well. 4 is recommened for isolated peak and 2 for slighted overlapped peak
        :param width_base: --optional, number of points to calculate the base counts, higher width_base reduces uncertainty but enlarges influence of irrelevant signals
        :param draw_peak_region: --optional, draw images of peak area

        :return: averagted total area and uncertainty
        '''

        list_peaks_distance=[abs(self.ergcal.chn2erg(peak.centroid)-erg) for peak in self.list_peaks]
        peak_index=list_peaks_distance.index(min(list_peaks_distance))
        peak=self.list_peaks[peak_index]

        peak_half_width=max(int(self.FWHM(self.ergcal.chn2erg(peak.centroid))/self.ergcal.c1/2.35*4),10)
        peak.peak_right_boundary-=3                                                                                                                                                                                            
        while True: # Determine the peak area boundary using minimum point and (s[i]-s[i+1])<s[i]**0.5
            count1=self.list_eff[peak.peak_right_boundary]
            count2=self.list_eff[peak.peak_right_boundary+1]
            count3=self.list_eff[peak.peak_right_boundary+2]
            count4=self.list_eff[peak.peak_right_boundary+3]
            # if (count2-count1)**2<=(count1) and (count3-count2)**2<=(count2) and (count4-count3)**2<=(count3):
            if (count2-count1)**2<=(count1) and (count3-count2)**2<=(count2):
                peak.peak_right_boundary+=1
                break # estimation of peak area boundary use counts difference (should be small than 1 sigma for succesive points)
            elif peak.peak_right_boundary<=peak.centroid+peak_half_width:
               peak.peak_right_boundary+=1
            else:
                minimum_right=sgn.argrelmin(np.array(self.list_eff[peak.centroid:min(peak.centroid+peak_half_width,len(self.list_eff))]))
                if len(minimum_right[0])!=0:
                    peak.peak_right_boundary=peak.centroid+minimum_right[0][0] # if not succed in limited window length, use rough estimation of peak area boundary with minumum points
                else:
                    peak.peak_right_boundary=peak.centroid+int(peak_half_width*0.75)
                break
        
        peak.peak_left_boundary+=3
        while True: # Determine the peak area boundary using minimum point and (s[i]-s[i+1])<s[i]**0.5
            count1=self.list_eff[peak.peak_left_boundary] 
            count2=self.list_eff[peak.peak_left_boundary-1]
            count3=self.list_eff[peak.peak_left_boundary-2]
            count4=self.list_eff[peak.peak_left_boundary-3]
            # if (count2-count1)**2<=(count1) and (count3-count2)**2<=(count2) and (count4-count3)**2<=(count3):
            if (count2-count1)**2<=(count1) and (count3-count2)**2<=(count2):
                peak.peak_left_boundary-=1
                break 
            elif peak.peak_left_boundary<=peak.centroid-peak_half_width:
                peak.peak_left_boundary-=1
            else:
                minimum_left=sgn.argrelmin(np.array(self.list_eff[max(peak.centroid-peak_half_width,0):peak.centroid+1]))
                if len(minimum_left[0])!=0:
                    peak.peak_left_boundary=peak.centroid-len(self.list_eff[max(peak.centroid-peak_half_width,0):peak.centroid+1])+1+minimum_left[0][-1]
                else:
                    peak.peak_left_boundary=peak.centroid-int(peak_half_width*0.75)
                break

        list_counts=[]
        list_variance=[]
        width_base_left=width_base
        width_base_right=width_base
        
        tolerance_base=1.5
        gap_width_left=0
        while True:
            all_background=1
            for chn in range(width_base_left+avg_number):
                count1=self.list_eff[peak.peak_left_boundary-gap_width_left-chn]
                count2=self.list_eff[peak.peak_left_boundary-gap_width_left-chn+1]
                if (count1-count2)**2<=tolerance_base*count1:
                    pass
                else:
                    all_background*=0
                    break
            if all_background==1:
                break
            else:
                gap_width_left+=1

        gap_width_right=0
        while True:
            all_background=1
            for chn in range(width_base_right+avg_number):
                count1=self.list_eff[min(peak.peak_right_boundary+gap_width_right+chn,len(self.list_eff)-1)]
                count2=self.list_eff[min(peak.peak_right_boundary+gap_width_right+chn+1,len(self.list_eff)-1)]
                if (count1-count2)**2<=tolerance_base*count1:
                    pass
                else:
                    all_background*=0
                    break
            if (all_background==1) or (peak.peak_right_boundary+gap_width_right+avg_number==len(self.list_eff)):
                break
            else:
                gap_width_right+=1

        # peak.peak_left_boundary=3824   

        for avg_i in range(avg_number): # average peak area of avg_number times area counting
            base_boundary_left=peak.peak_left_boundary-gap_width_left-avg_i-1
            base_boundary_right=peak.peak_right_boundary+gap_width_right+avg_i+1
            total_counts=sum(self.list_eff[peak.peak_left_boundary:peak.peak_right_boundary+1])
            
            if base_method==0: # linear background substraction
                base_counts_left=sum(self.list_eff[base_boundary_left-width_base_left+1:base_boundary_left+1])/width_base_left
                base_counts_right=sum(self.list_eff[base_boundary_right:base_boundary_right+width_base_right])/width_base_right
                base_total_counts=(base_counts_left+base_counts_right)/2*(peak.peak_right_boundary-peak.peak_left_boundary+1) # simulated base area using linear assumption
                variance_base=1/4*(peak.peak_right_boundary-peak.peak_left_boundary+1)**2*(base_counts_left/width_base_left+base_counts_right/width_base_right)
            elif base_method==1: # Polynomial curvefitting background substraction
                list_erg_short=self.list_erg[base_boundary_left-width_base_left:base_boundary_left]
                list_erg_short.extend(self.list_erg[base_boundary_right+1:base_boundary_right+width_base_right])
                list_eff_short=self.list_eff[base_boundary_left-width_base_left:base_boundary_left]
                list_eff_short.extend(self.list_eff[base_boundary_right+1:base_boundary_right+width_base_right])
                param,covar=np.polyfit(list_erg_short,list_erg_short,3,cov=True) # fitted parameters and covariance matrix
                param_covar=[covar[i][i] for i in range(len(param))]
                polyvals=np.polyval(param,self.list_erg)
                base_total_counts=sum(polyvals[peak.peak_left_boundary:peak.peak_right_boundary+1])
                variance_base=sum(np.polyval(param_covar,self.list_erg)[peak.peak_left_boundary:peak.peak_right_boundary+1])
            elif base_method==2: # SINP background substraction
                if not hasattr(self,'list_background'):
                    self.SNIP()
                base_total_counts=sum(self.list_background[peak.peak_left_boundary:peak.peak_right_boundary+1])
                variance_base=base_total_counts
            else:
                pass

            variance=variance_base+total_counts
            list_counts.append(total_counts-base_total_counts)
            list_variance.append(variance) # note the area and variance in each calculation
            # print(total_counts-base_total_counts,variance**0.5/(total_counts-base_total_counts)) # For refineing 
        net_counts_average=sum(list_counts)/avg_number
        variance_statistics=sum(list_variance)/avg_number # variance caused by counting statistics

        if avg_number!=1:
            variance_pamameter_selection=0
            for avg_i in range(avg_number):
                variance_pamameter_selection+=(list_counts[avg_i]-net_counts_average)**2
            variance_pamameter_selection=variance_pamameter_selection/(avg_number-1) # variance caused by parameter selection
        else:
            variance_pamameter_selection=0
        total_variance=variance_statistics
        # total_variance=variance_statistics+variance_pamameter_selection  # Generally not required
        uncertainty=total_variance**0.5/net_counts_average 
        base_start_left=peak.peak_left_boundary-gap_width_left
        base_start_right=peak.peak_right_boundary+gap_width_right
        return net_counts_average, uncertainty, total_counts, peak.peak_left_boundary, peak.peak_right_boundary,base_start_left,base_start_right

    def report(self):
        report='centroid energy net_counts_average, uncertainty, total_counts, peak.peak_left_boundary, peak.peak_right_boundary,base_start_left,base_start_right \n'
        for peak in self.list_peaks:
            net_counts_average, uncertainty, total_counts, peak_left_boundary, peak_right_boundary,base_start_left,base_start_right=self.areaTPA(self.ergcal.chn2erg(peak.centroid),width_base=3,base_method=2)
            report+='{:>8d} {:>8.2f} \t {:>10.2f} {:>8.2f}% \t {:>10.2f} {:>8d} {:>8d} {:>8d} {:>8d}\n'.format(peak.centroid,self.ergcal.chn2erg(peak.centroid),net_counts_average, uncertainty*100, total_counts, peak_left_boundary, peak_right_boundary,base_start_left,base_start_right)

        print(report)
        return 0
    def export_GV(self,filepath,sample_discription='No sample description was entered.',live_time=10000,total_time=10000):
    # def write(list_eff,filepath,increment=0,slope=0.29317,sample_discription='No sample description was entered.',live_time=10000,total_time=10000):
        '''
        Write spectrum files in Spe format. Requires a existing Spe file to dump file

        :param list_eff: --mandatory, counts per channel list of length 8192(for spe format)
        :param filepath: --mandatory, object filepath of spe file
        :param live_time: --optional, live measurement time, 10000 by default
        :param total time: --optional, total measurement time, 10000 by default
        :param sample_discription: --optional, description about sample in string format
        :return: spe formate spectrum in selected filepath
        '''
        spe_open=open(filepath,mode='w',encoding='utf-8')
        spe_open.write('$SPEC_ID:\n')
        spe_open.write(sample_discription+'\n')
        spe_open.write('$SPEC_REM: \n')
        spe_open.write('DET# 1 \n')
        spe_open.write('DETDESC# 1#GB \n')
        spe_open.write('AP# GammaVision Version 6.08 \n')
        spe_open.write('$DATE_MEA: \n')
        time_structure=time.localtime(time.time())
        spe_open.write('%02i/%02i/%04i %02i:%02i:%02i\n'%(time_structure[1],time_structure[2],time_structure[0],time_structure[3],time_structure[4],time_structure[5])) # Needs verification
        spe_open.write('$MEAS_TIM: \n')
        spe_open.write('%5.5s %5.5s\n'%(live_time,total_time))
        spe_open.write('$DATA: \n0 8191\n')    
        i=0
        while True:
            try:
                eff=int(self.list_eff_unsmoothed[i])
                spe_open.write('%8d\n'%(eff))
                i+=1
            except:
                spe_open.write('$ROI:\n')
                spe_open.write('0\n')
                spe_open.write('$PRESETS:\n')
                spe_open.write('None\n0\n0\n')
                spe_open.write('$ENER_FIT:\n')
                spe_open.write('%8.6f %8.6f\n'%(self.ergcal.c0,self.ergcal.c1))
                spe_open.write('$MCA_CAL:\n')
                spe_open.write('3\n')
                spe_open.write('%.6e %.6e %.6e keV\n'%(self.ergcal.c0,self.ergcal.c1,self.ergcal.c2))
                spe_open.write('$SHAPE_CAL:\n')
                spe_open.write('3\n')
                spe_open.write('6.418442E+000 0.000000E+000 0.000000E+000')
                spe_open.close()
                break
        return 0


if __name__ == '__main__':
    Cs137=Spectrum()
    Cs137.import_GV('E:\\Spectrum_analysis\\Spectrum_data\\GM20燃耗球\\2020-07-21-1800s-023-目标球（反）.chn')
    # Cs137.smooth_wavelet(level=10)
    Cs137.smooth_fourier(factor=0.5)
    Cs137.search_simple()
    Cs137.draw_spectrum()
    Cs137.report()