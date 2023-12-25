import os
import time
from attr import s
import numpy as np
import scipy.optimize as opt
import scipy.signal as sgn
import matplotlib.pyplot as plt

from math import ceil
from scipy.optimize import curve_fit

def EnergyCalibration(channel,slope=0.293676,increment=1.605054):
    '''
    Return corresponding energy for a given channel, with energy calibration slope and increment specified

    :param channel: --mandatory
    :param slope: --Optional
    :param increment --Optional
    :return energy in keV
    '''
    return channel*slope+increment

def read(filepath):
    '''
    Read spectrum files in Spe or MCNP output format.

    :param filepath: --filepath of spe or MCNP-output file
    :return: list_erg and list_eff
    '''
    list_erg=[] 
    list_eff=[]  
    with open(filepath,'r',encoding='utf-8') as fileopen:
        filelines=fileopen.readlines()
    indl=0 # index of file line
    try:
        if filepath[-4:] == '.spe':
            while True:
                line=filelines[indl]
                if '$ENER_FIT:' in line:
                    slope=float(filelines[indl+1].split()[1])
                    increment=float(filelines[indl+1].split()[0])
                    break
                else:
                    indl+=1
            indl=0
            while True:
                line=filelines[indl]
                indl+=1
                if '0 8191' in line:
                    break
            chn=0 # index of energy channel

            while True:
                try:
                    erg=EnergyCalibration(chn,slope=slope,increment=increment)
                    line=filelines[indl]
                    line=line.split()
                    list_eff.append(float(line[0]))
                    list_erg.append(erg)
                    chn+=1
                    indl+=1
                except:
                    del list_erg[-1]
                    del list_eff[-1]
                    break
        else:
            while True:
                line=filelines[indl]
                indl+=1
                if '1tally   8        ' in line:
                    indl+=6
                    break
            while True:
                try:
                    line=filelines[indl]
                    line=line.split()
                    list_erg.append(float(line[0])*1000)
                    list_eff.append(float(line[1]))
                    indl+=1
                except:
                    break
        return list_erg,list_eff
    except:
        print('Error-This file cannot be interpreted as Spe or MCNP-ouput format.')
        return 0


def write(list_eff,filepath,increment=0,slope=0.29317,sample_discription='No sample description was entered.',live_time=10000,total_time=10000,):
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
            eff=int(list_eff[i])
            spe_open.write('%8d\n'%(eff))
            i+=1
        except:
            spe_open.write('$ROI\n')
            spe_open.write('0\n')
            spe_open.write('$PRESETS\n')
            spe_open.write('None\n0\n0\n')
            spe_open.write('$ENER_FIT\n')
            spe_open.write('%8.6f %8.6f\n'%(increment,slope))
            spe_open.write('$MCA_CAL\n')
            spe_open.write('3\n')
            spe_open.write('0.000000E+000 0.000000E+000 0.000000E+000\n')
            spe_open.write('$SHAPE_CAL\n')
            spe_open.write('3\n')
            spe_open.write('0.000000E+000 0.000000E+000 0.000000E+000\n')
            spe_open.close()
            break


def SNIP(list_erg,list_eff):
    '''
    Background substration with SNIP method


    :param list_erg: --mandatory, energy list 
    :param list_eff: --mandatory, counts per channel list of same length with list_erg

    :return: spectrum of background in length of list_eff
    '''
    m=10
    l=len(list_erg)
    list_v=[0]*l
    list_background=[0]*l
    for i in range(l):
        list_v[i]=np.log(np.log((list_eff[i]+1)**0.5+1)+1)
    for p in range(1,m+1):
        for i in range(p,l-p):
            list_v[i]=min(list_v[i],0.5*(list_v[i-p]+list_v[i+p]))
    for i in range(l):
        list_background[i]=(np.exp(np.exp(list_v[i])-1))**2-1

    return list_background

def smooth(list_eff,width=21,order=4):
    '''
    Return smoothed spectrum with Savitzky-Golay filter

    :param list_eff: --mandatory, unsmoothed spectrum
    :param width: --Optional, window width, default 21, larger width leads to boardening
    :param order: --Optinal, Polynome order,default 4, higher order leads to stronger effect
    '''
    list_eff_smoothed=sgn.savgol_filter(list_eff,width,order)
    list_eff=list_eff_smoothed
    return list_eff

def areaTPA(list_erg,list_eff,peak,avg_number=4,width_base=4,draw_peak_region=0,base_method=0):
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
    list_eff_smoothed=[] # smoothing the spectrum- should find other options
    l=len(list_eff)
    for i in range(l):
        # list_eff_smoothed.append((2*list_eff[i]+list_eff[max(i-1,0)]+list_eff[min(i+1,l-1)])/4) # 3-points smoothing method  
        # list_eff_smoothed.append((4*list_eff[i]+list_eff[max(i-1,0)]+list_eff[min(i+1,l-1)])/6) # 3-points Simpson smoothing method 
        list_eff_smoothed.append((-3*list_eff[max(i-2,0)]+12*list_eff[max(i-1,0)]+17*list_eff[i]+12*list_eff[min(i+1,l-1)]-3*list_eff[min(i+2,l-1)])/35) # # 5-points smoothing method 
    list_eff=list_eff_smoothed

    width=25
    peak_boundary_right=peak
    while True: # Determine the peak area boundary using minimum point and (s[i]-s[i+1])<s[i]**0.5
        count1=list_eff[peak_boundary_right]
        count2=list_eff[peak_boundary_right+1]
        count3=list_eff[peak_boundary_right+2]
        count4=list_eff[peak_boundary_right+3]
        if (count2-count1)**2<=(count1) and (count3-count2)**2<=(count2) and (count4-count3)**2<=(count3):
            break # estimation of peak area boundary use counts difference (should be small than 1 sigma for succesive points)
        elif (peak_boundary_right-peak)<=width:
            peak_boundary_right+=1
        else:
            minimum_right=sgn.argrelmin(np.array(list_eff[peak:peak+min(width,len(list_eff)-peak)]))
            peak_boundary_right=peak+minimum_right[0][0] # if not succed in limited window length, use rough estimation of peak area boundary with minumum points
            break
    peak_boundary_left=peak
    while True: # Determine the peak area boundary using minimum point and (s[i]-s[i+1])<s[i]**0.5
        count1=list_eff[peak_boundary_left] 
        count2=list_eff[peak_boundary_left-1]
        count3=list_eff[peak_boundary_left-2]
        count4=list_eff[peak_boundary_left-3]
        if (count2-count1)**2<=(count1) and (count3-count2)**2<=(count2) and (count4-count3)**2<=(count3):
            break # estimation of peak area boundary use counts difference (should be small than 1 sigma for succesive points)
        elif (peak-peak_boundary_left)<=width:
            peak_boundary_left-=1
        else:
            minimum_left=sgn.argrelmin(np.array(list_eff[peak-min(width,peak):peak+1]))
            peak_boundary_left=peak+minimum_left[0][-1]-len(list_eff[peak-min(width,peak):peak+1])
            break # if not succed in limited window length, use rough estimation of peak area boundary with minumum points
    if draw_peak_region==1:
        filepath=os.path.abspath('.')
        if not os.path.exists(filepath+'\\areaTPA_peak_region'):os.mkdir(filepath+'\\areaTPA_peak_region')
        ymax=max(list_eff[peak_boundary_left:peak_boundary_right])
        plt.plot(list_eff[peak-width:peak+width+1],'r.')
        plt.axvline(x=peak_boundary_left-peak+width,ymax=list_eff[peak_boundary_left]/ymax,linestyle='--')
        plt.axvline(x=peak_boundary_right-peak+width,ymax=list_eff[peak_boundary_right]/ymax,linestyle='--')
        plt.xlabel('Channel')
        plt.ylabel('Counts')
        plt.title('Peak Region Plot\n Peak=%i Erg=%4.4s'%(peak,list_erg[peak]))

        plt.savefig(filepath+'\\areaTPA_peak_region\\peak=%i.jpg'%(peak))
        plt.close()


    list_counts=[]
    list_variance=[]
    width_base_left=width_base
    width_base_right=width_base
    for base_half_width in range(avg_number): # average peak area of avg_number times area counting
        base_boundary_left=peak_boundary_left-base_half_width-1
        base_boundary_right=peak_boundary_right+base_half_width+1
        total_counts=sum(list_eff[peak_boundary_left:peak_boundary_right+1])
        
        if base_method==0: # linear background substraction
            base_counts_left=sum(list_eff[base_boundary_left-width_base_left+1:base_boundary_left+1])/width_base_left
            base_counts_right=sum(list_eff[base_boundary_right:base_boundary_right+width_base_right])/width_base_right
            base_total_counts=(base_counts_left+base_counts_right)/2*(peak_boundary_right-peak_boundary_left+1) # simulated base area using linear assumption
            variance_base=1/4*(peak_boundary_right-peak_boundary_left+1)**2*(base_counts_left/width_base_left+base_counts_right/width_base_right)
        elif base_method==1: # Polynomial curvefitting background substraction
            list_erg_short=list_erg[base_boundary_left-width_base_left+1:base_boundary_left]
            list_erg_short.extend(list_erg[base_boundary_right:base_boundary_right+width_base_right])
            list_erg_short=list_eff[base_boundary_left-width_base_left+1:base_boundary_left]
            list_erg_short.extend(list_eff[base_boundary_right:base_boundary_right+width_base_right])
            param,covar=np.polyfit(list_erg_short,list_erg_short,3,cov=True) # fitted parameters and covariance matrix
            param_covar=[covar[i][i] for i in range(len(param))]
            polyvals=np.polyval(param,list_erg)
            base_total_counts=sum(polyvals[peak_boundary_left:peak_boundary_right+1])
            variance_base=sum(np.polyval(param_covar,list_erg)[peak_boundary_left:peak_boundary_right+1])
        elif base_method==2: # SINP background substraction
            list_base=SNIP(list_erg,list_eff)
            base_total_counts=sum(list_base[peak_boundary_left:peak_boundary_right+1])

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
        for base_half_width in range(avg_number):
            variance_pamameter_selection+=(list_counts[base_half_width]-net_counts_average)**2
        variance_pamameter_selection=variance_pamameter_selection/(avg_number-1) # variance caused by parameter selection
    else:
        variance_pamameter_selection=0
    total_variance=variance_statistics
    # total_variance=variance_statistics+variance_pamameter_selection  # Generally not required
    uncertainty=total_variance**0.5/net_counts_average 
    return net_counts_average, uncertainty,variance_statistics,variance_pamameter_selection


def FWHMfit(erg):
    '''
    This funcion provides characteristic HWHM value of GEM20 detector. 
    The fitting function is based on calibration spectrum.

    :param erg: --mandatory, energy of peak centroid in keV 

    :return: FWHM in keV 
    '''
    a=-0.00139;b=0.00517;c=-0.376;erg=erg/1000        
    return (a+b*pow(erg+c*erg**2,0.5))*1000

x=1
def areaCowell(list_erg,list_eff,peak,sigma=-1,N1=1,N2=3):
    '''
    This function realizes average Covell peak area counting method.Covell method has fixed peak area boundary and baseline boundary.

    '''
    if sigma==-1:
        sigma=FWHMfit(list_erg[peak])/2.355    
    Nmin=int(sigma*2.355/0.29317*N1)
    Nmax=int(sigma*2.355/0.29317*N2)
    list_cts_peak=[]
    list_err=[]
    for i in range(Nmin,Nmax+1):
        chn_l=peak-i;chn_r=peak+i
        cts_l=list_eff[chn_l]
        cts_r=list_eff[chn_r]
        cts_back=(cts_l+cts_r)/2*(chn_r-chn_l+1)
        cts_total=sum(list_eff[chn_l:chn_r+1])
        cts_peak=cts_total-cts_back
        list_cts_peak.append(cts_peak)
        err=pow(sum(list_eff[chn_l+1:chn_r])+(chn_r-chn_l-1)**2*(cts_r+cts_l)/4,0.5)
        list_err.append(err)
        # print('\nchn_l=%i,chn_r=%i,cts_back=%i,cts_peak=%i,err=%f,uncertainty=%f'%(chn_l,chn_r,cts_back,cts_peak,err,err/cts_peak))

        # a=(list_eff[peak-N]-list_eff[peak+N])/(-2*N)
        # b=list_eff[peak-N]-a*(peak-N)
        # list_back=[a*n+b for n in range(len(list_erg))]
        # plt.axvline(x=list_erg[peak],linestyle='--')
        # plt.axvline(x=list_erg[peak-N],linestyle='--')
        # plt.axvline(x=list_erg[peak+N],linestyle='--')
        # plt.plot(list_erg[peak-50:peak+50],list_eff[peak-50:peak+50])
        # plt.plot(list_erg[peak-50:peak+50],list_back[peak-50:peak+50])
        # plt.show()
    cts_peak_avg=sum(list_cts_peak)/(Nmax-Nmin+1)
    err_avg=sum(list_err)/(Nmax-Nmin+1)
    err2_avg=0
    for i in range(Nmax-Nmin+1):
        err2_avg+=(list_cts_peak[i]-cts_peak_avg)**2
    err2_avg=err2_avg/(Nmax-Nmin)
    err3_avg=pow(err_avg**2+err2_avg,0.5)
    return cts_peak_avg, err3_avg


def base(x,a,b,c):
    return a*x**2+b*x+c
def intbase(x1,x2,a,b,c):
    y1=a*x1**3/3+b*x1**2/2+c*x1
    y2=a*x2**3/3+b*x2**2/2+c*x2
    return y2-y1

def areaQuittner(list_erg,list_eff,peak,k_l,k_r,w_l,w_r,lmax):
    list_cts_peak=[]
    list_err=[]
    ind=0
    cts_total=sum(list_eff[peak-w_l:peak+w_l+1])
    for i in range(lmax):
        chn_l=peak-w_l-i-k_l
        chn_r=peak+w_r+i+k_l
        list_x=list_erg[chn_l-k_l:chn_l+k_l+1]
        list_x.extend(list_erg[chn_r-k_l:chn_r+k_l+1])
        list_y=list_eff[chn_l-k_l:chn_l+k_l+1]
        list_y.extend(list_eff[chn_r-k_l:chn_r+k_l+1])
        a,b,c=opt.curve_fit(base,list_x,list_y)[0]
        list_y2=[base(x,a,b,c) for x in list_x]
        # plt.plot(list_x,list_y)
        # plt.plot(list_x,list_y2)
        # plt.ylim([0,2000])
        # plt.show()
        x=opt.curve_fit(base,list_x,list_y)[1]
        x1=list_erg[peak-w_l]
        x2=list_erg[peak+w_l]
        cts_back=intbase(x1,x2,a,b,c)*(2*w_l+1)/(2*w_l)/0.29317
        cts_peak=cts_total-cts_back
        list_cts_peak.append(cts_peak)
        # err=pow(sum(list_eff[chn_l:chn_r+1])+(chn_r-chn_l+1)*(chn_r-chn_l)*(cts_r+cts_l)/4,0.5)
        err=0
        list_err.append(err)
        ind+=1
        # print('\nchn_l=%i,chn_r=%i,cts_back=%i,cts_peak=%i,err=%f,uncertainty=%f'%(chn_l,chn_r,cts_back,cts_peak,err,err/cts_peak))

        # a=(list_eff[peak-N]-list_eff[peak+N])/(-2*N)
        # b=list_eff[peak-N]-a*(peak-N)
        # list_back=[a*n+b for n in range(len(list_erg))]
        # plt.axvline(x=list_erg[peak],linestyle='--')
        # plt.axvline(x=list_erg[peak-N],linestyle='--')
        # plt.axvline(x=list_erg[peak+N],linestyle='--')
        # plt.plot(list_erg[peak-50:peak+50],list_eff[peak-50:peak+50])
        # plt.plot(list_erg[peak-50:peak+50],list_back[peak-50:peak+50])
        # plt.show()
    cts_peak_avg=sum(list_cts_peak)/(ind)
    err_avg=sum(list_err)/(ind)
    err2_avg=0
    for i in range(ind):
        err2_avg+=(list_cts_peak[i]-cts_peak_avg)**2
    err2_avg=err2_avg/(ind)
    err3_avg=pow(err_avg**2+err2_avg,0.5)
    return cts_peak_avg, err3_avg


def areaSterlinski(list_erg,list_eff,peak,sigma=-1,N=3,m=10):

    Nmin=int(sigma*2.355/0.29317*1)
    Nmax=int(sigma*2.355/0.29317*N)

    list_l=list_eff[peak-50:peak+1]
    list_r=list_eff[peak:peak+51]
    min_l=sgn.argrelmin(np.array(list_l))
    min_r=sgn.argrelmin(np.array(list_r))
    chn_l0=peak+min_l[0][-1]-50
    chn_r0=peak+min_r[0][0]

    list_cts_peak=[]
    list_err=[]
    ind=0
    for i in range(m):
        chn_l=chn_l0-i
        chn_r=chn_r0+i
        Nmax2=min(peak-chn_l,chn_r-peak)
        Nmax3=min(Nmax,Nmax2)
        for j in range(Nmin,Nmax3+1):
            chn_l2=peak-j;chn_r2=peak+j
            cts_l=list_eff[chn_l]
            cts_r=list_eff[chn_r]
            cts_back=(cts_l+cts_r)/2*(chn_r2-chn_l2+1)
            
            cts_total=sum(list_eff[chn_l2:chn_r2+1])
            cts_peak=cts_total-cts_back
            list_cts_peak.append(cts_peak)
            err=pow(sum(list_eff[chn_l:chn_r+1])+(chn_r-chn_l+1)*(chn_r-chn_l)*(cts_r+cts_l)/4,0.5)
            list_err.append(err)
            ind+=1
        # print('\nchn_l=%i,chn_r=%i,cts_back=%i,cts_peak=%i,err=%f,uncertainty=%f'%(chn_l,chn_r,cts_back,cts_peak,err,err/cts_peak))

        # a=(list_eff[peak-N]-list_eff[peak+N])/(-2*N)
        # b=list_eff[peak-N]-a*(peak-N)
        # list_back=[a*n+b for n in range(len(list_erg))]
        # plt.axvline(x=list_erg[peak],linestyle='--')
        # plt.axvline(x=list_erg[peak-N],linestyle='--')
        # plt.axvline(x=list_erg[peak+N],linestyle='--')
        # plt.plot(list_erg[peak-50:peak+50],list_eff[peak-50:peak+50])
        # plt.plot(list_erg[peak-50:peak+50],list_back[peak-50:peak+50])
        # plt.show()
    cts_peak_avg=sum(list_cts_peak)/(ind)
    err_avg=sum(list_err)/(ind)
    err2_avg=0
    for i in range(ind):
        err2_avg+=(list_cts_peak[i]-cts_peak_avg)**2
    err2_avg=err2_avg/(ind)
    err3_avg=pow(err_avg**2+err2_avg,0.5)
    return cts_peak_avg, err3_avg



# 半高宽拟合
def FWHMfit(erg):
    '''
    Function to fit the FWHM curve based on semi-empirical formula.

    :param erg: --mandatory, energy in keV
    :return FWHM in keV
    '''
    a=-0.00139;b=0.00517;c=-0.376;erg=erg/1000        
    return (a+b*pow(erg+c*erg**2,0.5))*1000

# 高斯函数拟合
def gauss(x,x0,sigma,h):
    y=h*np.exp(-(x-x0)**2/(2*sigma**2))
    return y
def baseline(x,a,b):
    return a*x+b
def gausspeak(erg,main_peak_energy,main_peak_uncertainty,main_peak_amplitude,baseline_slope,baseline_increment):
    return main_peak_amplitude*np.exp(-(erg-main_peak_energy)**2/(2*main_peak_uncertainty**2))+baseline_slope*erg+baseline_increment
def double_gausspeak(erg,main_peak_energy,main_peak_uncertainty,main_peak_amplitude,baseline_slope,baseline_increment,second_peak_energy,second_peak_uncertainty,second_peak_amplitude):
    y=main_peak_amplitude*np.exp(-(erg-main_peak_energy)**2/(2*main_peak_uncertainty**2))+baseline_slope*erg+baseline_increment
    y+=second_peak_amplitude*np.exp(-(erg-second_peak_energy)**2/(2*second_peak_uncertainty**2))
    return y
def triple_gausspeak(erg,main_peak_energy,main_peak_uncertainty,main_peak_amplitude,baseline_slope,baseline_increment,second_peak_energy,second_peak_uncertainty,second_peak_amplitude,third_peak_energy,third_peak_uncertainty,third_peak_amplitude):
    y=main_peak_amplitude*np.exp(-(erg-main_peak_energy)**2/(2*main_peak_uncertainty**2))+baseline_slope*erg+baseline_increment
    y+=second_peak_amplitude*np.exp(-(erg-second_peak_energy)**2/(2*second_peak_uncertainty**2))
    y+=third_peak_amplitude*np.exp(-(erg-third_peak_energy)**2/(2*third_peak_uncertainty**2))
    return y    

# 本底+高斯峰拟合
def peakfit(list_erg,list_eff,main_peak,window_width=20):
    list_erg_windowed=list_erg[max(main_peak-window_width,0):min(main_peak+window_width+1,len(list_erg))]
    list_eff_windowed=list_eff[max(main_peak-window_width,0):min(main_peak+window_width+1,len(list_eff))]
    main_peak_channel=main_peak-max(main_peak-window_width,0)
    main_peak_energy=list_erg_windowed[main_peak_channel]

    filepath=os.path.abspath('.')
    if not os.path.exists(filepath+'\\peakfit'):os.mkdir(filepath+'\\peakfit')
       
    list_erg_windowed_sorted_by_counts,list_eff_windowed_sorted_by_counts=zip(*sorted(zip(list_eff_windowed,list_erg_windowed))) # preliminary baseline estimation
    l=int(0.3*len(list_eff_windowed)) # assume the smallest 30% of signals are definitely baseline
    baseline_slope_0,baseline_increment_0=curve_fit(baseline,list_eff_windowed_sorted_by_counts[:l],list_erg_windowed_sorted_by_counts[:l])[0]
    main_peak_energy_0=list_erg_windowed[main_peak_channel] # estimation on main peak energy
    main_peak_uncertainty_0=abs(FWHMfit(main_peak_energy)/2.35) # estimation on main peak sigma broadening
    main_peak_amplitude_0=list_eff_windowed[main_peak_channel]-baseline_increment_0-baseline_slope_0*list_erg_windowed[main_peak_channel]
    # list_eff_windowed_prefitted=[gausspeak(erg,main_peak_energy_0,main_peak_uncertainty_0,main_peak_amplitude_0,baseline_slope_0,baseline_increment_0) for erg in list_erg_windowed]  

    try:
        para_set=curve_fit(triple_gausspeak,list_erg_windowed,list_eff_windowed,p0=[main_peak_energy_0,main_peak_uncertainty_0,main_peak_amplitude_0,baseline_slope_0,baseline_increment_0])[0]
        # main_peak_energy,main_peak_uncertainty,main_peak_amplitude,baseline_slope,baseline_increment=curve_fit(gausspeak,list_erg_windowed,list_eff_windowed,p0=[main_peak_energy_0,main_peak_uncertainty_0,main_peak_amplitude_0,baseline_slope_0,baseline_increment_0])[0]
        # list_eff_windowed_fitted=[gausspeak(erg,main_peak_energy,main_peak_uncertainty,main_peak_amplitude,baseline_slope,baseline_increment) for erg in list_erg_windowed]
        # plt.plot(list_erg_windowed,list_eff_windowed_fitted,label='Fitted')
        # plt.legend()
        # plt.savefig(filepath+'\\peakfit\\fitted_peak=%s_fitted.jpg'%(main_peak))
        # plt.close()
        return para_set
    except:

        return 'Unable to fit the peak at channel=%i'%(main_peak)
    
    

if __name__=="__main__":

    # list_erg,total_counts=read(r'C:\Users\Alber\AppData\Local\VirtualStore\Program Files (x86)\LANL\MCNP5\全能谱分析\Data'+'\\Cs-137.spe')
    # list_erg,total_counts=read(r'C:\Users\Alber\AppData\Local\VirtualStore\Program Files (x86)\LANL\MCNP5\全能谱分析\Data'+'\\Co-60.spe')
    # list_erg,total_counts=read(r'C:\Users\Alber\AppData\Local\VirtualStore\Program Files (x86)\LANL\MCNP5\全能谱分析\Data'+'\\Eu-152.spe')
    # list_erg,list_eff=read('E:\\MCNP_data\\Spectrum_analysis\\Data\\Eu-152.spe')
    list_erg,list_eff=read('E:\\MCNP_data\\Spectrum_analysis\\Data\\2020-07-21-1800s-023-目标球.spe')
    # list_peak_channel=[416,837,1176,2656,3288,3704,3792,4803]
    # # peakfit(list_erg,total_counts,837,window_width=40)
    
    # # report=open('report.txt','w')
    # # for peaks in list_peak_channel:
    # #     print(peaks)
    # #     report.write('peak=%i\t'%(peaks)+str(areaTPA(total_counts,peak=peaks,avg_number=3,width_base=6))+'\n')
    # # report.close
    # import numpy.fft as fourier
    # from scipy.signal import find_peaks
    # def fourier_smooth(list_eff,main_peak,window_width=20):
    #     list_erg_windowed=list_erg[max(main_peak-window_width,0):min(main_peak+window_width+1,len(list_erg))]
    #     list_eff_windowed=list_eff[max(main_peak-window_width,0):min(main_peak+window_width+1,len(list_eff))]
    #     main_peak_channel=main_peak-max(main_peak-window_width,0)
    #     main_peak_energy=list_erg_windowed[main_peak_channel]
    #     list_eff_transformed=fourier.fft(list_eff_windowed)

    #     limit_frequency=15;r=0.001
    #     list_eff_transformed_filtered_high_frequency=[]
    #     list_eff_transformed_filtered_low_frequency=[]
    #     for i in range(len(list_eff_transformed)):
    #         if (i <= limit_frequency-1) or (i >= len(list_eff_transformed)-limit_frequency-1):
    #             list_eff_transformed_filtered_low_frequency.append(list_eff_transformed[i])
    #             list_eff_transformed_filtered_high_frequency.append(list_eff_transformed[i]*r)
    #         else:
    #             list_eff_transformed_filtered_low_frequency.append(list_eff_transformed[i]*r)
    #             list_eff_transformed_filtered_high_frequency.append(list_eff_transformed[i])   
    #     list_eff_transformed_filtered_low_frequence_retransformed=fourier.ifft(list_eff_transformed_filtered_low_frequency)
    #     list_eff_transformed_filtered_high_frequence_retransformed=fourier.ifft(list_eff_transformed_filtered_high_frequency)
    #     return list_eff_transformed_filtered_low_frequence_retransformed,list_eff_transformed_filtered_high_frequence_retransformed

    # list_eff_transformed_filtered_low_frequence_retransformed,list_eff_transformed_filtered_high_frequence_retransformed=fourier_smooth(total_counts,list_peak_channel[1],window_width=40)
    # list_found_peaks=find_peaks(np.real(np.array(list_eff_transformed_filtered_low_frequence_retransformed)),prominence=50)


    # # print(len(list_found_peaks[0]))
    # # main_peak=list_peak_channel[1]
    # # window_width=40
    # # list_eff=total_counts[max(main_peak-window_width,0):min(main_peak+window_width+1,len(list_erg))]
    # # plt.plot(list_eff,label='Unfilterd')
    # # plt.plot(list_eff_transformed_filtered_high_frequence_retransformed,label='High_frequency')
    # # plt.plot(list_eff_transformed_filtered_low_frequence_retransformed,label='Low_frequency')
    # # for x in list_found_peaks[0]:
    # #     plt.axvline(x)
    # # plt.legend()
    # # plt.show()

    # print(len(list_found_peaks[0]))
    # main_peak=list_peak_channel[1]
    # window_width=40
    # list_erg_windowed=list_erg[max(main_peak-window_width,0):min(main_peak+window_width+1,len(list_erg))]
    # list_eff_windowed=np.real(list_eff_transformed_filtered_low_frequence_retransformed)
    # list_erg_windowed_sorted_by_counts,list_eff_windowed_sorted_by_counts=zip(*sorted(zip(list_eff_windowed,list_erg_windowed))) # preliminary baseline estimation
    # l=int(0.3*len(list_eff_windowed)) # assume the smallest 30% of signals are definitely baseline
    # baseline_slope_0,baseline_increment_0=curve_fit(baseline,list_eff_windowed_sorted_by_counts[:l],list_erg_windowed_sorted_by_counts[:l])[0]
        
    # p0=[[list_erg_windowed[list_found_peaks[0][3]],1,list_eff_windowed[list_found_peaks[0][3]]/0.293177,baseline_slope_0,baseline_increment_0,list_erg_windowed[list_found_peaks[0][1]],1,list_eff_windowed[list_found_peaks[0][1]]/0.293177,list_erg_windowed[list_found_peaks[0][2]],1,list_eff_windowed[list_found_peaks[0][2]]/0.293177]]
    # list_eff_pre_fitted=[triple_gausspeak(erg,list_erg_windowed[list_found_peaks[0][3]],1,list_eff_windowed[list_found_peaks[0][3]]-baseline_increment_0-baseline_slope_0*list_erg_windowed[list_found_peaks[0][3]],baseline_slope_0,baseline_increment_0,list_erg_windowed[list_found_peaks[0][1]],1,list_eff_windowed[list_found_peaks[0][1]]-baseline_increment_0-baseline_slope_0*list_erg_windowed[list_found_peaks[0][1]],list_erg_windowed[list_found_peaks[0][2]],1,list_eff_windowed[list_found_peaks[0][2]]-baseline_increment_0-baseline_slope_0*list_erg_windowed[list_found_peaks[0][2]]) for erg in list_erg_windowed]
    # para=curve_fit(triple_gausspeak,list_erg_windowed,list_eff_windowed,p0)[0]
    # list_eff_fitted=[triple_gausspeak(erg,para[0],para[1],para[2],para[3],para[4],para[5],para[6],para[7],para[8],para[9],para[10]) for erg in list_erg_windowed]
    # print(para)
    # list_eff=total_counts[max(main_peak-window_width,0):min(main_peak+window_width+1,len(list_erg))]
    # plt.plot(list_eff,label='Unfilterd')
    # plt.plot(list_eff_transformed_filtered_high_frequence_retransformed,label='High_frequency')
    # plt.plot(list_eff_transformed_filtered_low_frequence_retransformed,label='Low_frequency')
    # plt.plot(list_eff_pre_fitted,label='Pre_fitted')
    # plt.plot(list_eff_fitted,label='fitted')
    # for x in list_found_peaks[0]:
    #     plt.axvline(x)
    # plt.legend()
    # plt.show()



    # ind_peak=0
    # peak_num=len(list_found_peaks[0])
    # list_para_set=[[0,1,0,0,0]]*peak_num
    # para_set=peakfit(list_erg,list_eff,list_found_peaks[0][ind_peak],window_width=5)
    # list_para_set[ind_peak]=para_set
    # list_eff_new=[]
    # for i in range(len(list_eff)):
    #     eff_new=list_eff[i]
    #     for ind_para in range(peak_num):
    #         eff_new-=gausspeak(list_erg[i],list_para_set[ind_para][0],list_para_set[ind_para][1],list_para_set[ind_para][2],list_para_set[ind_para][3],list_para_set[ind_para][4])
    #     list_eff_new.append(eff_new)
    # plt.plot(list_eff,label='origin')
    # ind=1
    # plt.plot(list_eff_new,label='ind=%i'%(ind))




    # plt.legend()
    # plt.savefig(os.path.abspath('.'+'\\peak=%i-ind=%i.jpg'%(peak,ind)))
    # critia=0.1

    # while 1:
    #     ind+=1
    #     error=0
    #     for i in range(len(list_eff_new)):
    #         error+=(list_eff_new[i]/list_eff[i])**2
    #     if error<=critia: break
    #     ind_peak=(ind_peak+1)%peak_num
    #     para_set=peakfit(list_erg,list_eff_new,list_found_peaks[0][ind_peak],window_width=5)
    #     list_para_set[ind_peak]=para_set
    #     list_eff_new=[]
    #     for i in range(len(list_eff)):
    #         eff_new=list_eff[i]
    #         for ind_para in range(peak_num):
    #             eff_new-=gausspeak(list_erg[i],list_para_set[ind_para][0],list_para_set[ind_para][1],list_para_set[ind_para][2],list_para_set[ind_para][3],list_para_set[ind_para][4])
    #         list_eff_new.append(eff_new)
    #     plt.plot(list_eff_new,label='ind=%i'%(ind))
    #     plt.legend()
    #     plt.savefig(os.path.abspath('.'+'\\peak=%i-ind=%i.jpg'%(peak,ind)))


    # plt.plot(list_eff_transformed_filtered_high_frequency,'x',label='High_frequency')
    # plt.plot(list_eff_transformed_filtered_low_frequency,'x',label='Low_frequency')
    # plt.plot(list_eff_transformed_abs,'x',label='Total')
    # plt.yscale('log')



