import numpy as np
import matplotlib.pyplot as plt
import func_smooth as sm
import func_area as area

def gauss(centroid, stderror, area, channels):
    channels_peak = np.zeros(shape=channels.size)
    channels_peak[max(int(centroid-4*stderror), 0):min(int(centroid+4*stderror), channels.shape[0])] = \
        channels[max(int(centroid-4*stderror), 0):min(int(centroid+4*stderror), channels.shape[0])].copy()
    amplitude = area / stderror / (2*np.pi)**0.5
    return np.exp( -(channels-centroid)**2 / stderror**2 /2 ) * amplitude 

def SNIP(list_eff, m=25):
    
    list_eff = np.array(list_eff)
    loglog_eff = np.log(np.log((list_eff+1)**0.5+1)+1)
    loglog_eff = np.pad(loglog_eff, pad_width=m, mode='edge')

    for p in range(m+1, 1, -1):
        for i in range(p,loglog_eff.shape[0]-p):
            loglog_eff[i]=min(loglog_eff[i],(loglog_eff[i-p]+loglog_eff[i+p])/2)
    baseline = (np.exp(np.exp(loglog_eff)-1)-1)**2-1
    baseline = baseline[m: -m]
    return baseline

def fit(list_eff, peak_info):
    baseline = SNIP(list_eff)
    net_spectrum = list_eff - baseline
    area = []
    # plt.plot(spectrum, '--', label='original')
    # plt.plot(net_spectrum, '--', label='SNIP')
    # plt.legend()
    # plt.show()
    for peak in peak_info:
        peak_centroid, peak_sigma, _ = peak
        gross = sum(net_spectrum[int(peak_centroid-3*peak_sigma): int(peak_centroid+3*peak_sigma)])
        noise_left = sum(net_spectrum[int(peak_centroid-4*peak_sigma): int(peak_centroid-3*peak_sigma)]) / len(net_spectrum[int(peak_centroid-4*peak_sigma): int(peak_centroid-3*peak_sigma)])
        noise_right = sum(net_spectrum[int(peak_centroid+3*peak_sigma): int(peak_centroid+4*peak_sigma)]) / len(net_spectrum[int(peak_centroid+3*peak_sigma): int(peak_centroid+4*peak_sigma)])
        area.append( gross - (noise_left+noise_right)/2 * len(net_spectrum[int(peak_centroid-3*peak_sigma): int(peak_centroid+3*peak_sigma)]) )
    return net_spectrum
  
channels = np.arange(0, 200, step=1)
peak_info = [[30, 2.5, 1000], [70, 3, 500],[110, 4,200], 
             [150, 4.5, 600], [165, 4.3, 300]] 
peaks = []
baseline = (channels/300 - 1)**2 * 5 + 20
noise = baseline**0.3 * np.random.normal(size=channels.shape)
spectrum = baseline + noise
for i, item in enumerate(peak_info):
    peak = gauss(peak_info[i][0], peak_info[i][1], peak_info[i][2], channels)
    peaks.append(peak)
    spectrum += peak

centroid5 = sm.centroid(spectrum, L=2)
savitzy52 = sm.savol(spectrum, L=2, R=2)
fourier_smoothed_low = sm.fourier(spectrum, f=0.3, sigma=2, type='low')
fourier_smoothed_gauss = sm.fourier(spectrum, f=0.3, sigma=1, type='gauss')
fourier_smoothed_mixed = sm.fourier(spectrum, f=0.3, sigma=1, type='mixed')

spectrums = [spectrum, centroid5, savitzy52, fourier_smoothed_low, fourier_smoothed_gauss, fourier_smoothed_mixed]
spectrums_names = ['original spectrum', 'Centroid (N=5)', 'Savitzy-Golay (N=5, R=2)', 'Fourier (ω=0.3)', 'Gauss-Fourier (σ=1)', 'Mixed-Fourier( ω=0.3, σ=1)']

# for s in spectrums:
#     print(fit(s, peak_info))
                   
def smoothingfactor(unsmoothed, smoothed, N=9):
    numerator = np.zeros(shape=unsmoothed[N:].shape)
    denominator = np.zeros(shape=unsmoothed[N:].shape)
    for i in range(N):
        numerator += smoothed[i:len(smoothed)-N+i]**2
        denominator += unsmoothed[i:len(smoothed)-N+i]**2
    return numerator /  denominator

N=9
# for i,item in enumerate(spectrums):
#     plt.plot(item+np.full(shape=item.shape, fill_value=i*10), label=spectrums_names[i], marker='.', linestyle='--', markersize=3, linewidth=1)
#     # plt.plot(smoothingfactor(spectrum, item, N)+np.full(shape=item[N:].shape, fill_value=i), label=spectrums_names[i], marker='.', linestyle='--', markersize=3, linewidth=1)
#     SNR = sum(item[30-10:30+10]**2) / sum(item[30+10:30+30]**2)
#     print(SNR)
# plt.xlabel('Channels',fontdict={'fontsize':10,'family':'times new roman','style':'normal'})
# plt.ylabel('Counts',fontdict={'fontsize':10,'family':'times new roman','style':'normal'})
# plt.title('(a) Comparaison of smoothing algorithms',fontdict={'fontsize':10,'family':'times new roman','style':'normal'})
# xticks=np.arange(0,200, 10)
# xlabels=['{:5.0f}'.format(tick) for tick in xticks]
# plt.xticks(xticks,xlabels,fontdict={'fontsize':8,'family':'times new roman','style':'normal','rotation':0})
# yticks=np.arange(0, 250, 20)
# ylabels=['{:5.1f}'.format(tick) for tick in yticks]
# plt.yticks(yticks,ylabels,fontdict={'fontsize':8,'family':'times new roman','style':'normal','rotation':0})
# plt.grid(which='both',axis='both',color='grey',linewidth=0.2,linestyle='-.')
# plt.legend()
# plt.legend(loc='upper right',fontsize='xx-small',prop={'family':'times new roman'},frameon=1,fancybox=0)
# plt.show()



# import numpy as np     
# def Ieff_savgol(L, R):
#     A = np.array( [ [ i**j for j in range(R+1) ] for i in range(-L, L+1) ] )
#     M = np.matmul(A, np.matmul( np.linalg.inv(np.matmul(A.T, A)), A.T ) )[L]
#     return sum(M**2)
# R=5
# for L in range(int((R+1)/2),5):
#     print(L, R, Ieff_savgol(L,R))
    
# import numpy as np  
# import matplotlib.pyplot as plt 
# def Ieff_fourier(w, s, type):
#     frequence = np.linspace(0, 100, 100)
#     gauss_transformed = lambda w: np.exp( -s**2 * w**2 / 2)
#     low_transformed = lambda i: 1 - (i>(100*w)) * (i<(100-100*w))
#     def mixed_transformed(spectrum):
#         for i, v in enumerate(spectrum):
#             if (i>(100*w)) and (i<(100-100*w)):
#                 j = -abs(i-50) + 50 - 100*w 
#                 spectrum[i] = np.exp( -s**2 * j**2 / 2)
#             else:
#                 spectrum[i] = 1
#         return spectrum
#     if type == 'gauss':
#         transformed = gauss_transformed(frequence)
#     elif type == 'low':
#         transformed = low_transformed(frequence) 
#     else:
#         transformed = mixed_transformed(frequence)
#     untransformed = np.fft.ifft(transformed)
#     plt.plot(untransformed, label=str(w)+str(s))
#     Ieff = sum(np.real(untransformed)**2)
#     return Ieff
# # for s in range(1,9):
# #     print(s, Ieff_fourier(0.1, s, type='gauss'))
# # for w in np.linspace(0.1, 0.4, 10):
# #     print(w, Ieff_fourier(w, 0.1, type='low')) 
# for w in np.linspace(0.1, 0.4, 2):
#     for s in range(1, 5):
#         print(w, s, Ieff_fourier(w, s, type='mixed'))  
# plt.legend()
# plt.show()

import scipy.signal as sgn
import numpy as np


def areaATPA(spectrum, erg,avg_number=4,width_base=3,base_method=0):



    list_peaks_distance=[abs(instance_spectrum.ergcal.chn2erg(peak.centroid)-erg) for peak in instance_spectrum.list_peaks]
    peak_index=list_peaks_distance.index(min(list_peaks_distance))
    peak=instance_spectrum.list_peaks[peak_index]

    peak_half_width=max(int(instance_spectrum.fwhmcal.FWHM(instance_spectrum.ergcal.chn2erg(peak.centroid))/instance_spectrum.ergcal.c1/2.35*4),10)
    peak.peak_right_boundary-=3                                                                                                                                                                                            
    while True: # Determine the peak area boundary using minimum point and (s[i]-s[i+1])<s[i]**0.5
        count1=instance_spectrum.list_eff[peak.peak_right_boundary]
        count2=instance_spectrum.list_eff[peak.peak_right_boundary+1]
        count3=instance_spectrum.list_eff[peak.peak_right_boundary+2]
        # count4=instance_spectrum.list_eff[peak.peak_right_boundary+3]
        # if (count2-count1)**2<=(count1) and (count3-count2)**2<=(count2) and (count4-count3)**2<=(count3):
        if (count2-count1)**2<=(count1) and (count3-count2)**2<=(count2):
            peak.peak_right_boundary+=1
            break # estimation of peak area boundary use counts difference (should be small than 1 sigma for succesive points)
        elif peak.peak_right_boundary<=peak.centroid+peak_half_width:
            peak.peak_right_boundary+=1
        else:
            minimum_right=sgn.argrelmin(np.array(instance_spectrum.list_eff[peak.centroid:min(peak.centroid+peak_half_width,len(instance_spectrum.list_eff))]))
            if len(minimum_right[0])!=0:
                peak.peak_right_boundary=peak.centroid+minimum_right[0][0] # if not succed in limited window length, use rough estimation of peak area boundary with minumum points
            else:
                peak.peak_right_boundary=peak.centroid+int(peak_half_width*0.75)
            break
    
    peak.peak_left_boundary+=3
    while True: # Determine the peak area boundary using minimum point and (s[i]-s[i+1])<s[i]**0.5
        count1=instance_spectrum.list_eff[peak.peak_left_boundary] 
        count2=instance_spectrum.list_eff[peak.peak_left_boundary-1]
        count3=instance_spectrum.list_eff[peak.peak_left_boundary-2]
        # count4=instance_spectrum.list_eff[peak.peak_left_boundary-3]
        # if (count2-count1)**2<=(count1) and (count3-count2)**2<=(count2) and (count4-count3)**2<=(count3):
        if (count2-count1)**2<=(count1) and (count3-count2)**2<=(count2):
            peak.peak_left_boundary-=1
            break 
        elif peak.peak_left_boundary<=peak.centroid-peak_half_width:
            peak.peak_left_boundary-=1
        else:
            minimum_left=sgn.argrelmin(np.array(instance_spectrum.list_eff[max(peak.centroid-peak_half_width,0):peak.centroid+1]))
            if len(minimum_left[0])!=0:
                peak.peak_left_boundary=peak.centroid-len(instance_spectrum.list_eff[max(peak.centroid-peak_half_width,0):peak.centroid+1])+1+minimum_left[0][-1]
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
            count1=instance_spectrum.list_eff[peak.peak_left_boundary-gap_width_left-chn]
            count2=instance_spectrum.list_eff[peak.peak_left_boundary-gap_width_left-chn+1]
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
            count1=instance_spectrum.list_eff[min(peak.peak_right_boundary+gap_width_right+chn,len(instance_spectrum.list_eff)-1)]
            count2=instance_spectrum.list_eff[min(peak.peak_right_boundary+gap_width_right+chn+1,len(instance_spectrum.list_eff)-1)]
            if (count1-count2)**2<=tolerance_base*count1:
                pass
            else:
                all_background*=0
                break
        if (all_background==1) or (peak.peak_right_boundary+gap_width_right+avg_number==len(instance_spectrum.list_eff)):
            break
        else:
            gap_width_right+=1

    # peak.peak_left_boundary=3824   

    for avg_i in range(avg_number): # average peak area of avg_number times area counting
        base_boundary_left=peak.peak_left_boundary-gap_width_left-avg_i-1
        base_boundary_right=peak.peak_right_boundary+gap_width_right+avg_i+1
        total_counts=sum(instance_spectrum.list_eff[peak.peak_left_boundary:peak.peak_right_boundary+1])
        
        if base_method==0: # linear background substraction
            base_counts_left=sum(instance_spectrum.list_eff[base_boundary_left-width_base_left+1:base_boundary_left+1])/width_base_left
            base_counts_right=sum(instance_spectrum.list_eff[base_boundary_right:base_boundary_right+width_base_right])/width_base_right
            base_total_counts=(base_counts_left+base_counts_right)/2*(peak.peak_right_boundary-peak.peak_left_boundary+1) # simulated base area using linear assumption
            variance_base=1/4*(peak.peak_right_boundary-peak.peak_left_boundary+1)**2*(base_counts_left/width_base_left+base_counts_right/width_base_right)
        elif base_method==1: # Polynomial curvefitting background substraction
            list_erg_short=instance_spectrum.list_erg[base_boundary_left-width_base_left:base_boundary_left]
            list(list_erg_short).extend(instance_spectrum.list_erg[base_boundary_right+1:base_boundary_right+width_base_right])
            list_eff_short=instance_spectrum.list_eff[base_boundary_left-width_base_left:base_boundary_left]
            list(list_eff_short).extend(instance_spectrum.list_eff[base_boundary_right+1:base_boundary_right+width_base_right])
            param,covar=np.polyfit(list_erg_short,list_erg_short,3,cov=True) # fitted parameters and covariance matrix
            param_covar=[covar[i][i] for i in range(len(param))]
            polyvals=np.polyval(param,instance_spectrum.list_erg)
            base_total_counts=sum(polyvals[peak.peak_left_boundary:peak.peak_right_boundary+1])
            variance_base=sum(np.polyval(param_covar,instance_spectrum.list_erg)[peak.peak_left_boundary:peak.peak_right_boundary+1])
        elif base_method==2: # SINP background substraction
            if not hasattr(instance_spectrum,'list_background'):
                instance_spectrum.SNIP()
            base_total_counts=sum(instance_spectrum.list_background[peak.peak_left_boundary:peak.peak_right_boundary+1])
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
    if net_counts_average!=0:
        uncertainty=total_variance**0.5/net_counts_average
    else:
        uncertainty=0
    base_start_left=peak.peak_left_boundary-gap_width_left
    base_start_right=peak.peak_right_boundary+gap_width_right
    
    import matplotlib.pyplot as plt
    list_IF=[(instance_spectrum.list_eff[i]-instance_spectrum.list_eff[i+1])**2/(instance_spectrum.list_eff[i]+1) for i in range(len(instance_spectrum.list_eff)-1)]
    
    # fig=plt.subplots(2,1)[1]
    # ax=fig[0]
    # ax.plot(instance_spectrum.list_erg[base_start_left:base_start_right],instance_spectrum.list_eff[base_start_left:base_start_right])
    # ax.plot(instance_spectrum.list_erg[peak.peak_left_boundary:peak.peak_right_boundary],instance_spectrum.list_eff[peak.peak_left_boundary:peak.peak_right_boundary])
    # ax.set_xilm(instance_spectrum.list_erg[peak.peak_left_boundary-10],instance_spectrum.list_erg[peak.peak_right_boundary+10])
    # ax=fig[1]
    # ax.plot(instance_spectrum.list_erg[peak.peak_left_boundary-10:peak.peak_right_boundary+10],list_IF[peak.peak_left_boundary-10:peak.peak_right_boundary+10],'x')
    # plt.show()

    return net_counts_average, uncertainty, total_counts, peak.peak_left_boundary, peak.peak_right_boundary,base_start_left,base_start_right

plt.plot(spectrum, label='original spectrum', marker='.', linestyle='--', markersize=3, linewidth=1)
plt.plot(fit(centroid5, peak_info), label='SNIP strippd', marker='.', linestyle='--', markersize=3, linewidth=1)
plt.xlabel('Channels',fontdict={'fontsize':10,'family':'times new roman','style':'normal'})
plt.ylabel('Counts',fontdict={'fontsize':10,'family':'times new roman','style':'normal'})
plt.title('(b) SNIP stripped spectrum',fontdict={'fontsize':10,'family':'times new roman','style':'normal'})
xticks=np.arange(0,200, 10)
xlabels=['{:5.0f}'.format(tick) for tick in xticks]
plt.xticks(xticks,xlabels,fontdict={'fontsize':8,'family':'times new roman','style':'normal','rotation':0})
yticks=np.arange(0, 200, 20)
ylabels=['{:5.1f}'.format(tick) for tick in yticks]
plt.yticks(yticks,ylabels,fontdict={'fontsize':8,'family':'times new roman','style':'normal','rotation':0})
plt.grid(which='both',axis='both',color='grey',linewidth=0.2,linestyle='-.')
plt.legend()
plt.legend(loc='upper right',fontsize='xx-small',prop={'family':'times new roman'},frameon=1,fancybox=0)
plt.show()