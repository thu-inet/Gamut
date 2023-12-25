import numpy as np
import os
from numpy import polyfit
import sko
import matplotlib.pyplot as plt


def fit1_background(list_erg,list_eff):
    zip_list=zip(list_eff,list_erg)
    zip_list_sort=sorted(zip_list)
    list_eff,list_erg=zip(*zip_list_sort)
    param=polyfit(list_erg[:4],list_eff[:4],1,cov=False)
    def fit1(x):
        return param[0]*x+param[1]
    return fit1

def fit2_peak(list_erg,list_eff):
    list_lneff=[np.log(max(y,1)) for y in list_eff]
    param,var=polyfit(list_erg,list_lneff,2,cov=True)
    sigma=pow(-1/2/param[0],0.5)
    erg0=-param[1]/param[0]/2
    peak=np.exp(param[2]+erg0**2/2/sigma**2)
    def fit2(x):
        return peak*np.exp(-(x-erg0)**2/2/sigma**2)
    return fit2

list_eff=[38,28,42,45,56,80,89,105,124,135,162,148,133,107,83,62,53,40,37,33,30]
list_erg=np.arange(0,len(list_eff),1)
fit1=fit1_background(list_erg,list_eff)
list_eff=[y-fit1(x) for x,y in zip(list_erg,list_eff)]
fit2=fit2_peak(list_erg,list_eff)
list_fit=[fit1(x)+fit2(x) for x in list_erg]
print(list_fit)

# # 高斯函数拟合
# def gauss(x,x0,sigma,h):
#     y=h*np.exp(-(x-x0)**2/(2*sigma**2))
#     return y
# def gausspeak(erg,main_peak_energy,main_peak_uncertainty,main_peak_amplitude,background_slope,background_increment):
#     return main_peak_amplitude*np.exp(-(erg-main_peak_energy)**2/(2*main_peak_uncertainty**2))+background_slope*erg+background_increment
# def double_gausspeak(erg,main_peak_energy,main_peak_uncertainty,main_peak_amplitude,background_slope,background_increment,second_peak_energy,second_peak_uncertainty,second_peak_amplitude):
#     y=main_peak_amplitude*np.exp(-(erg-main_peak_energy)**2/(2*main_peak_uncertainty**2))+background_slope*erg+background_increment
#     y+=second_peak_amplitude*np.exp(-(erg-second_peak_energy)**2/(2*second_peak_uncertainty**2))
#     return y
# def triple_gausspeak(erg,main_peak_energy,main_peak_uncertainty,main_peak_amplitude,background_slope,background_increment,second_peak_energy,second_peak_uncertainty,second_peak_amplitude,third_peak_energy,third_peak_uncertainty,third_peak_amplitude):
#     y=main_peak_amplitude*np.exp(-(erg-main_peak_energy)**2/(2*main_peak_uncertainty**2))+background_slope*erg+background_increment
#     y+=second_peak_amplitude*np.exp(-(erg-second_peak_energy)**2/(2*second_peak_uncertainty**2))
#     y+=third_peak_amplitude*np.exp(-(erg-third_peak_energy)**2/(2*third_peak_uncertainty**2))
#     return y    

# # 本底+高斯峰拟合
# def peakfit(list_erg,list_eff,main_peak,window_width=20):
#     # list_erg_windowed=list_erg[max(main_peak-window_width,0):min(main_peak+window_width+1,len(list_erg))]
#     # list_eff_windowed=list_eff[max(main_peak-window_width,0):min(main_peak+window_width+1,len(list_eff))]
#     main_peak_channel=main_peak-max(main_peak-window_width,0)
#     main_peak_energy=list_erg_windowed[main_peak_channel]

#     filepath=os.path.abspath('.')
#     if not os.path.exists(filepath+'\\peakfit'):os.mkdir(filepath+'\\peakfit')
       
#     list_erg_windowed_sorted_by_counts,list_eff_windowed_sorted_by_counts=zip(*sorted(zip(list_eff_windowed,list_erg_windowed))) # preliminary background estimation
#     l=int(0.3*len(list_eff_windowed)) # assume the smallest 30% of signals are definitely background
#     background_slope_0,background_increment_0=curve_fit(background,list_eff_windowed_sorted_by_counts[:l],list_erg_windowed_sorted_by_counts[:l])[0]
#     main_peak_energy_0=list_erg_windowed[main_peak_channel] # estimation on main peak energy
#     main_peak_uncertainty_0=abs(FWHMfit(main_peak_energy)/2.35) # estimation on main peak sigma broadening
#     main_peak_amplitude_0=list_eff_windowed[main_peak_channel]-background_increment_0-background_slope_0*list_erg_windowed[main_peak_channel]
#     # list_eff_windowed_prefitted=[gausspeak(erg,main_peak_energy_0,main_peak_uncertainty_0,main_peak_amplitude_0,background_slope_0,background_increment_0) for erg in list_erg_windowed]  

#     try:
#         para_set=curve_fit(triple_gausspeak,list_erg_windowed,list_eff_windowed,p0=[main_peak_energy_0,main_peak_uncertainty_0,main_peak_amplitude_0,background_slope_0,background_increment_0])[0]
#         # main_peak_energy,main_peak_uncertainty,main_peak_amplitude,background_slope,background_increment=curve_fit(gausspeak,list_erg_windowed,list_eff_windowed,p0=[main_peak_energy_0,main_peak_uncertainty_0,main_peak_amplitude_0,background_slope_0,background_increment_0])[0]
#         # list_eff_windowed_fitted=[gausspeak(erg,main_peak_energy,main_peak_uncertainty,main_peak_amplitude,background_slope,background_increment) for erg in list_erg_windowed]
#         # plt.plot(list_erg_windowed,list_eff_windowed_fitted,label='Fitted')
#         # plt.legend()
#         # plt.savefig(filepath+'\\peakfit\\fitted_peak=%s_fitted.jpg'%(main_peak))
#         # plt.close()
#         return para_set
#     except:

#         return 'Unable to fit the peak at channel=%i'%(main_peak)

# def background_P1(chn,slope,increment):
#     return slope*chn+increment

# def gausspeak(erg,main_peak_energy,main_peak_uncertainty,main_peak_amplitude):
#     y=main_peak_amplitude*np.exp(-(erg-main_peak_energy)**2/(2*main_peak_uncertainty**2))
#     return y

# # 半高宽拟合
# def FWHMfit(erg):
#     '''
#     Function to fit the FWHM curve based on semi-empirical formula.

#     :param erg: --mandatory, energy in keV
#     :return FWHM in keV
#     '''
#     a=-0.00139;b=0.00517;c=-0.376;erg=erg/1000        
#     return (a+b*pow(erg+c*erg**2,0.5))*1000    

# # 本底+高斯峰拟合
# def peakfit(list_erg,list_eff,main_peak,chn_left,chn_right):
#     list_erg_windowed=list_erg[chn_left:chn_right]
#     list_eff_windowed=list_eff[chn_left:chn_right]
#     main_peak_channel=main_peak-max(main_peak-chn_left,0)

       
#     list_erg_windowed_sorted_by_counts,list_eff_windowed_sorted_by_counts=zip(*sorted(zip(list_eff_windowed,list_erg_windowed))) # preliminary background estimation
#     l=int(0.2*(chn_right-chn_left)) # assume the smallest 30% of signals are definitely background
#     background_slope_0,background_increment_0=curve_fit(background_P1,list_eff_windowed_sorted_by_counts[:l],list_erg_windowed_sorted_by_counts[:l])[0]

#     main_peak_energy_0=list_erg_windowed[main_peak_channel] # estimation on main peak energy
#     main_peak_uncertainty_0=abs(FWHMfit(main_peak_energy_0)/2.35) # estimation on main peak sigma broadening
#     main_peak_amplitude_0=list_eff_windowed[main_peak_channel]-background_increment_0-background_slope_0*list_erg_windowed[main_peak_channel]
#     # list_eff_windowed_prefitted=[gausspeak(erg,main_peak_energy_0,main_peak_uncertainty_0,main_peak_amplitude_0,background_slope_0,background_increment_0) for erg in list_erg_windowed]  

#     try:
#         para_set=curve_fit(triple_gausspeak,list_erg_windowed,list_eff_windowed,p0=[main_peak_energy_0,main_peak_uncertainty_0,main_peak_amplitude_0,background_slope_0,background_increment_0])[0]
#         # main_peak_energy,main_peak_uncertainty,main_peak_amplitude,background_slope,background_increment=curve_fit(gausspeak,list_erg_windowed,list_eff_windowed,p0=[main_peak_energy_0,main_peak_uncertainty_0,main_peak_amplitude_0,background_slope_0,background_increment_0])[0]
#         # list_eff_windowed_fitted=[gausspeak(erg,main_peak_energy,main_peak_uncertainty,main_peak_amplitude,background_slope,background_increment) for erg in list_erg_windowed]
#         # plt.plot(list_erg_windowed,list_eff_windowed_fitted,label='Fitted')
#         # plt.legend()
#         # plt.savefig(filepath+'\\peakfit\\fitted_peak=%s_fitted.jpg'%(main_peak))
#         # plt.close()
#         return para_set
#     except:

        # return 'Unable to fit the peak at channel=%i'%(main_peak)

# filepath=r'D:\User\Spectra\2020-07-21-1800s-023-目标球.Spe'

# list_erg,list_eff=read.read(filepath)
# chn_left=2716
# chn_right=2750
# list_erg_windowed=list_erg[chn_left:chn_right]
# list_eff_windowed=list_eff[chn_left:chn_right]
# plt.plot(list_eff_windowed,label='origin')

# list_erg_windowed_sorted_by_counts,list_eff_windowed_sorted_by_counts=zip(*sorted(zip(list_erg_windowed,list_eff_windowed))) # preliminary background estimation
# l=int(0.2*(chn_right-chn_left)) # assume the smallest 30% of signals are definitely background
# background_slope_0,background_increment_0=curve_fit(background_P1,list_erg_windowed_sorted_by_counts[:l],list_eff_windowed_sorted_by_counts[:l])[0]
# list_background=[]
# for i in range(len(list_erg_windowed)):
#     list_background.append(background_P1(list_erg_windowed[i],background_slope_0,background_increment_0))
#     list_eff_windowed[i]=list_eff_windowed[i]-background_P1(list_erg_windowed[i],background_slope_0,background_increment_0)
# # plt.plot(list_background,label='background')
# # plt.plot(list_eff_windowed,label='no background')

# para_set_1=curve_fit(gausspeak,list_erg_windowed,list_eff_windowed,p0=[801.87,1,100],maxfev=10000)[0]
# print(para_set_1)
# list_eff_windowed_fitted_1=[]
# main_peak_energy,main_peak_uncertainty,main_peak_amplitude=para_set_1
# for i in range(len(list_erg_windowed)):
#     list_eff_windowed_fitted_1.append(gausspeak(list_erg_windowed[i],main_peak_energy,main_peak_uncertainty,main_peak_amplitude))
#     list_eff_windowed[i]=list_eff_windowed[i]-gausspeak(list_erg_windowed[i],main_peak_energy,main_peak_uncertainty,main_peak_amplitude)
# plt.plot(list_eff_windowed_fitted_1,label='fit1')
# plt.plot(list_eff_windowed,label='after fit1')

# para_set_2=curve_fit(gausspeak,list_erg_windowed,list_eff_windowed,p0=[800,1,100],maxfev=1000)[0]
# list_eff_windowed_fitted_2=[]
# main_peak_energy,main_peak_uncertainty,main_peak_amplitude=para_set_2
# # for i in range(len(list_erg_windowed)):
# #     list_eff_windowed[i]-=gausspeak(list_erg_windowed[i],main_peak_energy,main_peak_uncertainty,main_peak_amplitude)
# for i in range(len(list_erg_windowed)):
#     list_eff_windowed_fitted_2.append(gausspeak(list_erg_windowed[i],main_peak_energy,main_peak_uncertainty,main_peak_amplitude))
# plt.plot(list_eff_windowed_fitted_2,label='fit2')
# plt.legend()
# plt.show()