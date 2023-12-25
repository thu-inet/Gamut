import matplotlib.pyplot as plt
import SNIP
import scipy.signal as sgn
import os
import numpy as np



def areaTPA(list_erg,list_eff,peak,avg_number=4,width_base=4,base_method=0):
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

    width=40
    LEFF=len(list_eff)
    peak_boundary_right=peak+5
    while True: # Determine the peak area boundary using minimum point and (s[i]-s[i+1])<s[i]**0.5
        count1=list_eff[peak_boundary_right]
        count2=list_eff[peak_boundary_right+1]
        count3=list_eff[peak_boundary_right+2]
        count4=list_eff[peak_boundary_right+3]
        if (count2-count1)**2<=(count1) and (count3-count2)**2<=(count2) and (count4-count3)**2<=(count3):
            peak_boundary_right+=1
            break # estimation of peak area boundary use counts difference (should be small than 1 sigma for succesive points)
        elif (peak_boundary_right-peak)<=width:
            peak_boundary_right+=1
        else:
            minimum_right=sgn.argrelmin(np.array(list_eff[peak:peak+min(width,LEFF-peak)]))
            peak_boundary_right=peak+minimum_right[0][0] # if not succed in limited window length, use rough estimation of peak area boundary with minumum points
            break
    
    peak_boundary_left=peak-5
    while True: # Determine the peak area boundary using minimum point and (s[i]-s[i+1])<s[i]**0.5
        count1=list_eff[peak_boundary_left] 
        count2=list_eff[peak_boundary_left-1]
        count3=list_eff[peak_boundary_left-2]
        count4=list_eff[peak_boundary_left-3]
        if (count2-count1)**2<=(count1) and (count3-count2)**2<=(count2) and (count4-count3)**2<=(count3):
            peak_boundary_left-=1
            break # estimation of peak area boundary use counts difference (should be small than 1 sigma for succesive points)
        elif (peak-peak_boundary_left)<=width:
            peak_boundary_left-=1

        else:
            minimum_left=sgn.argrelmin(np.array(list_eff[peak-min(width,peak):peak+1]))
            peak_boundary_left=peak+minimum_left[0][-1]-len(list_eff[peak-min(width,peak):peak+1])
            break # if not succed in limited window length, use rough estimation of peak area boundary with minumum points

    list_counts=[]
    list_variance=[]
    width_base_left=width_base
    width_base_right=width_base
    
    tolerance_base=1.5
    gap_width_left=0
    while True:
        all_background=1
        for chn in range(width_base_left+avg_number):
            if (list_eff[peak_boundary_left-gap_width_left-chn]-list_eff[peak_boundary_left-gap_width_left-chn+1])**2<=tolerance_base*list_eff[peak_boundary_left-gap_width_left-chn]:
                all_background*=1
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
            if (list_eff[peak_boundary_right+gap_width_right+chn]-list_eff[peak_boundary_right+gap_width_right+chn-1])**2<=tolerance_base*list_eff[peak_boundary_right+gap_width_right+chn]:
                pass
            else:
                all_background*=0
                break
        if all_background==1:
            break
        else:
            gap_width_right+=1
    # gap_width_left=50
    # gap_width_right=60 
    for avg_i in range(avg_number): # average peak area of avg_number times area counting
        base_boundary_left=peak_boundary_left-gap_width_left-avg_i-1
        base_boundary_right=peak_boundary_right+gap_width_right+avg_i+1
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
        for avg_i in range(avg_number):
            variance_pamameter_selection+=(list_counts[avg_i]-net_counts_average)**2
        variance_pamameter_selection=variance_pamameter_selection/(avg_number-1) # variance caused by parameter selection
    else:
        variance_pamameter_selection=0
    total_variance=variance_statistics
    # total_variance=variance_statistics+variance_pamameter_selection  # Generally not required
    uncertainty=total_variance**0.5/net_counts_average 


    # filepath=os.path.abspath('.')
    # if not os.path.exists(filepath+'\\areaTPA_peak_region'):os.mkdir(filepath+'\\areaTPA_peak_region')
    # ymax=max(list_eff[peak_boundary_left:peak_boundary_right])
    # plt.plot(list_erg[peak-width:peak+width+1],list_eff[peak-width:peak+width+1],'b-',label='Smoothed')
    # plt.plot(list_erg[peak_boundary_left:peak_boundary_right+1],list_eff[peak_boundary_left:peak_boundary_right+1],'r-',label='Peak ROI')
    # plt.axvline(x=list_erg[peak_boundary_left],ymax=list_eff[peak_boundary_left]/ymax,linestyle='--')
    # plt.axvline(x=list_erg[peak_boundary_right],ymax=list_eff[peak_boundary_right]/ymax,linestyle='--')
    # plt.xlabel('Channel')
    # plt.ylabel('Counts')
    # plt.title('Peak Region Plot\n Peak=%i Erg=%4.4s'%(peak,list_erg[peak]))
    # plt.legend()
    # plt.savefig(filepath+'\\areaTPA_peak_region\\peak=%i.jpg'%(peak))
    # plt.close()

    return net_counts_average, uncertainty, total_counts, peak_boundary_left, peak_boundary_right, peak_boundary_left-gap_width_left,peak_boundary_right+gap_width_right

# list_eff_smoothed=[] # smoothing the spectrum
#     l=len(list_eff)
#     for i in range(l):
#         # list_eff_smoothed.append((2*list_eff[i]+list_eff[max(i-1,0)]+list_eff[min(i+1,l-1)])/4) # 3-points smoothing method  
#         # list_eff_smoothed.append((4*list_eff[i]+list_eff[max(i-1,0)]+list_eff[min(i+1,l-1)])/6) # 3-points Simpson smoothing method 
#         list_eff_smoothed.append((-3*list_eff[max(i-2,0)]+12*list_eff[max(i-1,0)]+17*list_eff[i]+12*list_eff[min(i+1,l-1)]-3*list_eff[min(i+2,l-1)])/35) # # 5-points smoothing method 
    #     width=40
    # if draw_peak_region==1:
    #     plt.plot(list_erg[peak-width:peak+width+1],list_eff[peak-width:peak+width+1],'g.',label='Unsmoothed')
    # list_eff=list_eff_smoothed
#     if draw_peak_region==1:
#         filepath=os.path.abspath('.')
#         if not os.path.exists(filepath+'\\areaTPA_peak_region'):os.mkdir(filepath+'\\areaTPA_peak_region')
#         ymax=max(list_eff[peak_boundary_left:peak_boundary_right])
#         plt.plot(list_erg[peak-width:peak+width+1],list_eff[peak-width:peak+width+1],'b-',label='Smoothed')
#         plt.plot(list_erg[peak_boundary_left:peak_boundary_right+1],list_eff[peak_boundary_left:peak_boundary_right+1],'r-',label='Peak ROI')
#         plt.axvline(x=list_erg[peak_boundary_left],ymax=list_eff[peak_boundary_left]/ymax,linestyle='--')
#         plt.axvline(x=list_erg[peak_boundary_right],ymax=list_eff[peak_boundary_right]/ymax,linestyle='--')
#         plt.xlabel('Channel')
#         plt.ylabel('Counts')
#         plt.title('Peak Region Plot\n Peak=%i Erg=%4.4s'%(peak,list_erg[peak]))
#         plt.legend()
#         plt.savefig(filepath+'\\areaTPA_peak_region\\peak=%i.jpg'%(peak))
#         plt.close()
