import scipy.signal as sgn
import numpy as np
def areaTPA(list_eff,peak,m=10):
    list_l=list_eff[peak-50:peak+1]
    list_r=list_eff[peak:peak+51]
    min_l=sgn.argrelmin(np.array(list_l))
    min_r=sgn.argrelmin(np.array(list_r))
    try:
        chn_l0=peak+min_l[0][-1]-50
    except:
        chn_l0=peak-50
    chn_r0=peak+min_r[0][0]
    list_cts_peak=[]
    list_err=[]
    for i in range(m):
        chn_l=chn_l0-i;chn_r=chn_r0+i
        cts_l=list_eff[chn_l]
        cts_r=list_eff[chn_r]
        cts_back=(cts_l+cts_r)/2*(chn_r-chn_l+1)
        cts_total=sum(list_eff[chn_l:chn_r+1])
        cts_peak=cts_total-cts_back
        list_cts_peak.append(cts_peak)
        err=pow(sum(list_eff[chn_l+1:chn_r])+(chn_r-chn_l-1)**2*(cts_r+cts_l)/4,0.5)
        list_err.append(err)
        # print('\nchn_l=%i,chn_r=%i,cts_back=%i,cts_peak=%i,err=%f,uncertainty=%f'%(chn_l,chn_r,cts_back,cts_peak,err,err/cts_peak))
        
        # plt.plot(list_eff[peak-50:peak+50])
        # plt.axvline(x=chn_l+50-peak)
        # plt.axvline(x=chn_r+50-peak)
        # plt.show()
    cts_peak_avg=sum(list_cts_peak)/(m)
    err_avg=sum(list_err)/(m)
    err2_avg=0
    for i in range(m):
        err2_avg+=(list_cts_peak[i]-cts_peak_avg)**2
    err2_avg=err2_avg/(m-1)
    err3_avg=pow(err_avg**2+err2_avg,0.5)
    return cts_peak_avg, err3_avg
