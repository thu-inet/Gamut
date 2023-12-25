'''
Smoothing
-definition: Smoothing is considered as an essentail primary operation in spectrum precessing. 
It uses synthetical infomation (e.g. counting data from multiple channels or the full spectrum) 
to reconstruct the counting data of each channel to compensate the statistical error in the spectrum. 

-function: Smoothing renders the spetrum smoother so that peaks can be easily detected and the peak area more acurately calculted.

-methods: There are several methods to do spectrum smoothing, including Savatizky-Golay, Fourier transformation, etc.
'''
def R(list_eff,list_eff_smthed):
    s1=0
    s2=0
    for i in len(1,list_eff-1):
        s1+=(list_eff[i]-list_eff[i-1])**2
        s2+=(list_eff_smthed[i]-list_eff_smthed[i-1])**2
    return s2/s1

def rho(list_eff,list_eff_smthed):
    mean1=sum(list_eff)/len(list_eff)
    mean2=sum(list_eff_smthed)/len(list_eff_smthed)
    rho=0
    for i in range(len(list_eff)):
        rho+=(list_eff[i]-mean1)*(list_eff_smthed[i]-mean2)/list_eff[i]**0.5/list_eff_smthed[i]**0.5
    return rho/(len(list_eff)-1)

def RMSE(list_eff,list_eff_smthed):
    rmse=0
    for i in range(len(list_eff)):
        rmse+=(list_eff_smthed[i]-list_eff[i])**2
    return (rmse/len(list_eff))**0.5

import numpy as np
def Per(list_eff,list_eff_smthed):
    mat_eff=np.mat(list_eff).T*np.mat(list_eff)
    engval,engvec=np.linalg.eigh(mat_eff)
    mat_eff_smthed=np.mat(list_eff_smthed).T*np.mat(list_eff_smthed)
    engval_smthed,engvec_smthed=np.linalg.eigh(mat_eff_smthed)
    return max(engval_smthed)/max(engval)

def Mean(list_eff,list_eff_smthed,mc,Mc):
    N001=sum(list_eff[mc:Mc+1])
    N002=sum(list_eff_smthed[mc:Mc+1])
    mean1=0
    mean2=0
    for i in range(-2,3):
        for j in range(-2,3):
            mean1+=abs(sum(list_eff[mc+i:Mc+1+j])-N001)
    for i in range(-2,3):
        for j in range(-2,3):
            mean2+=abs(sum(list_eff_smthed[mc+i:Mc+1+j])-N002)
    mean1*=1/24/N001
    mean2*=1/24/N002
    return mean2/mean1

def Std(list_eff,list_eff_smthed,mc,Mc):
    N001=sum(list_eff[mc:Mc+1])
    N002=sum(list_eff_smthed[mc:Mc+1])
    list_N1=[]
    list_N2=[]
    for i in range(-2,3):
        for j in range(-2,3):
            list_N1.append(abs(sum(list_eff[mc+i:Mc+1+j])-N001))
    for i in range(-2,3):
        for j in range(-2,3):
            list_N2.append(abs(sum(list_eff_smthed[mc+i:Mc+1+j])-N002))
    mean_N1=sum(list_N1)/len(list_N1)
    mean_N2=sum(list_N2)/len(list_N2)
    std1=0
    for N in list_N1:
        std1+=(N-mean_N1)**2
    std2=0
    for N in list_N2:
        std2+=(N-mean_N2)**2
    return std2/std1

def eval(list_eff,list_eff_smthed):
    Rx=R(list_eff,list_eff_smthed)
    rhox=rho(list_eff,list_eff_smthed)
    RMSEx=RMSE(list_eff,list_eff_smthed)
    Perx=Per(list_eff,list_eff_smthed)
    return Rx,rhox,RMSEx,Perx

def centroid(list_eff,winlen=3):
    list_eff_smthed=[eff for eff in list_eff]
    if winlen==3:
        for i in range(1,len(list_eff)-1):
            list_eff_smthed[i]=(list_eff[i-1]+2*list_eff[i]+list_eff[i+1])/4
    elif winlen==5:
        for i in range(2,len(list_eff)-2):
            list_eff_smthed[i]=(list_eff[i-2]+4*list_eff[i-1]+6*list_eff[i]+4*list_eff[i+1]+list_eff[i+2])/16
    elif winlen==7:
        for i in range(3,len(list_eff)-3):
            list_eff_smthed[i]=(list_eff[i-3]+6*list_eff[i-2]+15*list_eff[i-1]+20*list_eff[i]+15*list_eff[i+1]+6*list_eff[i+2]+list_eff[i+3])/64
    else:
        for i in range(1,len(list_eff)-1):
            list_eff_smthed[i]=(list_eff[i-1]+2*list_eff[i]+list_eff[i+1])/4
    return list_eff_smthed


from scipy.fftpack import fft, ifft
def fourier(list_eff,type='low',factor=0.5):
    trsed=fft(list_eff)
    trsed_smthed=trsed
    if type=='low':
        l=int(len(list_eff)/2*(1-factor))
        for i in range(l,len(trsed_smthed)-l):
            trsed_smthed[i]=0
    elif type=='gauss':
        factor=factor/len(trsed_smthed)/2
        for i in range(len(trsed_smthed)):
           trsed_smthed[i]*=np.exp(-factor*i**2)
    smthed=ifft(trsed_smthed)
    return smthed

from scipy.signal import savgol_filter
def savgol(list_eff,width,order):
    '''
    Simple Savgol
    '''
    smthed=savgol_filter(list_eff,width,order)
    return smthed

import numpy as np
def savgol2(list_eff,width,order):
    '''
    Svagol
    Ab+E=x
    -b=[b0,b1,b2,...,bn].n=order of polymomial
    -A=[[1,m,m^2,m^3,...,m^n],[1,m-1,(m-1)^2,)(m-1)^3,...,(m-1)^n],...,[1,-m,(-m)^2,-(m)^3,...,(-m)^n]].m=half width of transformation window
    -x=[y(m),y(m-1),...,y(0),...,y(-m)]

    solve: b=(At*A)^(-1)*At*x
    x=Ab=A*(At*A)^(-1)*At*x

    B=A*(At*A)^(-1)*At, x(new)=B*x=B*(x1,x2,...,xm)t
    xi(new)=B(ith row)*x

    '''
    half=int((width-1)/2)
    coefA=[]
    for i in range(-half,half+1):
        row=[i**j for j in range(order+1)]
        coefA.append(row)
    coefA=np.mat(coefA)
    coefB=(coefA*(coefA.T*coefA).I)*coefA.T
    coefB_x0=coefB[half]
    smthed=[]
    for i in range(half):
        smthed.append(list_eff[0])
    for i in range(half,len(list_eff)-half):
        window=coefB_x0*np.mat(list_eff[i-half:i+half+1]).T
        smthed.append(window[0,0])
    for i in range(half):
        smthed.append(list_eff[-1])
    return smthed

def convol(list_eff,type='gauss',winlen=7,funlen=4):
    list_eff_smthed=[eff for eff in list_eff]
    halfwin=int(winlen/2+1)
    if type=='gauss':
        list_convol=[np.exp(-4*np.ln(2)*j**2/funlen**2) for j in range(-halfwin,halfwin+1)]
        for i in range(halfwin,len(list_eff)-halfwin):
            y=0
            for j in (-halfwin,halfwin+1):
                y+=list_eff[i+j]*list_convol[j]
            list_eff_smthed[i]=y
    elif type=='cauchy':
        list_convol=[funlen**2/(funlen**2+4*j**2) for j in range(-halfwin,halfwin+1)]
        for i in range(halfwin,len(list_eff)-halfwin):
            y=0
            for j in (-halfwin,halfwin+1):
                y+=list_eff[i+j]*list_convol[j]
            list_eff_smthed[i]=y    
    elif type=='cos2':
        list_convol=[(np.cos(0.5*np.pi*j/funlen))**2 for j in range(-halfwin,halfwin+1)]
        for i in range(halfwin,len(list_eff)-halfwin):
            y=0
            for j in (-halfwin,halfwin+1):
                y+=list_eff[i+j]*list_convol[j]
            list_eff_smthed[i]=y
    elif type=='sech':
        list_convol=[2/(np.exp(2.634*j/funlen)+np.exp(-2.634*j/funlen)) for j in range(-halfwin,halfwin+1)]
        for i in range(halfwin,len(list_eff)-halfwin):
            y=0
            for j in (-halfwin,halfwin+1):
                y+=list_eff[i+j]*list_convol[j]
            list_eff_smthed[i]=y
    else:
        list_convol=[np.exp(-4*np.ln(2)*j**2/funlen**2) for j in range(-halfwin,halfwin+1)]
        for i in range(halfwin,len(list_eff)-halfwin):
            y=0
            for j in (-halfwin,halfwin+1):
                y+=list_eff[i+j]*list_convol[j]
            list_eff_smthed[i]=y
    return list_eff_smthed
    

    
if __name__=='__main__':
    import os
    import old.read as read
    import smoothing
    import numpy as np
    import matplotlib.pyplot as plt

    filepath=os.path.abspath('.')
    list_erg,list_eff=read.read(filepath+'\\Test_sample_spectrum.spe')
    mc=2000
    Mc=2200
    list_erg=list_erg[mc:Mc]
    list_eff=list_eff[mc:Mc]




    plt.plot(list_erg,list_eff,label='Unsmoothed')
    for factor in np.linspace(0.5,1,5):
        list_eff_smthed=smoothing.fourier(list_eff,factor)
        plt.plot(list_erg,list_eff_smthed,linewidth=0.3,label='Fourier factor=%s'%(factor))
    plt.yscale('log')
    plt.legend()
    plt.savefig(filepath+'\\Smoothing_Fourier.jpg')
    plt.close()


    # plt.plot(list_erg,list_eff,label='Unsmoothed')
    # for order in range(3,5):
    #     for width in range(5,10,2):
    #         # list_eff_smthed=smoothing.savgol(list_eff,width,order)
    #         # plt.plot(list_erg,list_eff_smthed,label='Savgol W=%i O=%i'%(width,order))
    #         list_eff_smthed=smoothing.savgol2(list_eff,width,order)
    #         plt.plot(list_erg,list_eff_smthed,linewidth=0.3,label='Savgol2 W=%i O=%i'%(width,order))
    # plt.yscale('log')
    # plt.legend()
    # plt.savefig(filepath+'\\Smoothing_Savgol.jpg')
    # plt.close()


    list_eff_smthed=smoothing.savgol2(list_eff,5,3)
        