import scipy.signal as sgn
import numpy as np


def areaTPA(instance_spectrum,erg,avg_number=4,width_base=3,base_method=0):

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