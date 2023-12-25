import spectrumfit
import scipy.optimize as opt
import numpy as np
import matplotlib.pyplot as plt
import scipy.signal as sgn

def areaSterlinski(list_erg,list_eff,peak,sigma=-1,N=3,m=10):
    if sigma ==-1:
        sigma=spectrumfit.spectrumfit(list_erg,list_eff,peak)['sigma']
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

def areaCowell(list_erg,list_eff,peak,sigma=-1,N1=1,N2=3):
    if sigma ==-1:
        sigma=spectrumfit.spectrumfit(list_erg,list_eff,peak)['sigma']
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