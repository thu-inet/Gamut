import spectrumfit
import scipy.signal as sgn
import numpy as np

def areaWasson(list_erg,list_eff,peak,sigma=-1,N=3,m=10):
    if sigma ==-1:
        sigma=spectrumfit.spectrumfit(list_erg,list_eff,peak)['sigma']
    Nmin=int(sigma*2.355/0.29317*1)
    Nmax=int(sigma*2.355/0.29317*N)
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
    ind=0
    for i in range(m):
        chn_l=chn_l0-i
        chn_r=chn_r0+i
        cts_l=list_eff[chn_l]
        cts_r=list_eff[chn_r]
        Nmax2=min(peak-chn_l,chn_r-peak)
        Nmax3=min(Nmax,Nmax2)
        for j in range(Nmin,Nmax3+1):
            chn_l2=peak-j;chn_r2=peak+j
            cts_back=(cts_l*(chn_r-peak)+cts_r*(peak-chn_l))*(2*j+1)/(chn_r-chn_l)
            cts_total=sum(list_eff[chn_l2:chn_r2+1])
            cts_peak=cts_total-cts_back
            list_cts_peak.append(cts_peak)
            err=pow(sum(list_eff[chn_l2:chn_r2+1])+(cts_l*(chn_r-peak)**2+cts_r*(peak-chn_l)**2)*(2*j+1)**2/(chn_r-chn_l)**2,0.5)
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

def areaWasson1(list_erg,list_eff,peak,n,m=10):
    list_l=list_eff[peak-50:peak+1]
    list_r=list_eff[peak:peak+51]
    min_l=sgn.argrelmin(np.array(list_l))
    min_r=sgn.argrelmin(np.array(list_r))
    try:
        chn_l0=peak+min_l[0][-1]-50
    except:
        chn_l0=peak-50
    chn_r0=peak+min_r[0][0]
    chn_l2=peak-n;chn_r2=peak+n
    cts_total=sum(list_eff[chn_l2:chn_r2+1])
    list_cts_peak=[]
    list_err=[]
    ind=0
    for i in range(m):
        chn_l=chn_l0-i
        chn_r=chn_r0+i
        cts_l=list_eff[chn_l]
        cts_r=list_eff[chn_r]
        cts_back=(cts_l*(chn_r-peak)+cts_r*(peak-chn_l))*(2*n+1)/(chn_r-chn_l)
        cts_peak=cts_total-cts_back
        list_cts_peak.append(cts_peak)
        err=pow(sum(list_eff[chn_l2:chn_r2+1])+(cts_l*(chn_r-peak)**2+cts_r*(peak-chn_l)**2)*(2*n+1)**2/(chn_r-chn_l)**2,0.5)
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

def areaWasson2(list_erg,list_eff,peak,sigma=-1,N=2.5,base_l=10,base_r=10):
    if sigma ==-1:
        sigma=spectrumfit.spectrumfit(list_erg,list_eff,peak)['sigma']
    Nmin=int(sigma*2.355/0.29317*1)
    Nmax=int(sigma*2.355/0.29317*N)
    chn_l=peak-base_l
    chn_r=peak+base_r    
    Nmax2=min(base_l,base_r)
    Nmax3=min(Nmax,Nmax2)
    list_cts_peak=[]
    list_err=[]
    ind=0
    for j in range(Nmin,Nmax3+1):
        chn_l2=peak-j;chn_r2=peak+j
        cts_l=list_eff[chn_l]
        cts_r=list_eff[chn_r]
        cts_back=(cts_l*(chn_r-peak)+cts_r*(peak-chn_l))*(2*j+1)/(chn_r-chn_l)
        cts_total=sum(list_eff[chn_l2:chn_r2+1])
        cts_peak=cts_total-cts_back
        list_cts_peak.append(cts_peak)
        err=pow(sum(list_eff[chn_l2:chn_r2+1])+(cts_l*(chn_r-peak)**2+cts_r*(peak-chn_l)**2)*(2*j+1)**2/(chn_r-chn_l)**2,0.5)
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