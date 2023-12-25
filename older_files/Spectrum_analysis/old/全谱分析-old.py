# 导入模块和设置
import matplotlib.pyplot as plt
import json
import numpy as np
import scipy.signal as sgn
import scipy.optimize as opt
import os
# filename=input('Please enter the file name=')
filename='Co-60.spe'
foldpath=r'C:\Users\Alber\AppData\Local\VirtualStore\Program Files (x86)\LANL\MCNP5\全能谱分析\\'
filepath=foldpath+filename
libpath=r'C:\Users\Alber\AppData\Local\VirtualStore\Program Files (x86)\LANL\MCNP5\Lib0.json'

# 读数据
def spectrumread(filepath):
    list_erg=[] 
    list_eff=[]  
    with open(filepath,'r',encoding='utf-8') as fileopen:
        filelines=fileopen.readlines()
    index=0
    if filepath[-4:] == '.spe':
        while True:
            line=filelines[index]
            index+=1
            if '0 8191' in line:
                break
        chn=0
        while True:
            try:
                erg=0.29317*chn
                line=filelines[index]
                line=line.split()
                list_eff.append(float(line[0]))
                list_erg.append(erg)
                chn+=1
                index+=1
            except:
                del list_erg[-1]
                del list_eff[-1]
                break
    else:
        while True:
            line=filelines[index]
            index+=1
            if '1tally   8        ' in line:
                index+=6
                break
        while True:
            try:
                line=filelines[index]
                line=line.split()
                list_erg.append(float(line[0])*1000)
                list_eff.append(float(line[1]))
                index+=1
            except:
                break
    return list_erg,list_eff
# 读核数据库
def loadlib(libpath):
    with open(libpath,'r',encoding='utf-8') as libopen:
        lib=json.load(libopen)
    return lib
# 画图
def spectrumplot(axe,list_erg,list_eff,title='',xlogscale=0,ylogscale=0,para='k-',linewidth=0.2):
    axe.plot(list_erg,list_eff,para,linewidth)
    if (ylogscale!=0):
        axe.set_yscale('log')
    if (xlogscale!=0):
        axe.set_xscale('log')
    axe.set_title('Spectrum'+title,fontsize=14)
    axe.set_xlabel('Energy (keV)',fontsize=12)
    axe.set_ylabel('Counts (cps) or FP Efficiency (1)',fontsize=12)
    axe.tick_params(axis='both',which='major',labelsize=10)
    # axe.set_xlim(100,2000)
    try:
        maxeff=max(list_eff)
        ymax=pow(10,int(np.log10(maxeff)+1))
        mineff=min(list_eff)
        ymin=pow(10,int(np.log10(mineff+1e-7)-1))
        axe.set_ylim(ymin,ymax)
    except:
        print('Some error happend in the ymax detemination')
    
# 寻峰
def findpeaks(list_erg,list_eff,prominence=200,width=5,distance=50,peaklocateerror=0.005):
    array_eff=np.array(list_eff)
    peaks,poperties=sgn.find_peaks(array_eff,distance,width,prominence)
    list_peaks=[]
    plt.plot(list_erg,list_eff)
    for peak in peaks:
        list_locatedpeaks=[]

        plt.plot(list_erg[peak-5:peak+5],list_eff[peak-5:peak+5],'r-')
        for libpeak in lib:
            error=abs((list_erg[peak]/libpeak['Energy']-1))
            if error<peaklocateerror:
                list_locatedpeaks.append(libpeak) 
        list_peaks.append(list_locatedpeaks)   
        # print('\nEnergy=',list_erg[peak],'keV,\tCandidate peaks=',len(list_peaks[i]))
        # for k in range(len(list_peaks[i])):
        #     print('\tPeak=',list_peaks[i][k])
        
        # plt.plot(list_erg[peak-50:peak+50],list_eff[peak-50:peak+50])
        # for k in range(len(list_peaks[i])):
        #     plt.axvline(x=list_peaks[i][k]['Energy'],linestyle='--')
        # plt.savefig(foldpath+r'\find_peaks\peak='+str(peak)[:5]+'_Energy='+str(peak*0.29317)[:5]+'.jpg')
        # plt.close()

    plt.savefig(foldpath+'\\find_peaks\\'+'peaknumber=%s_prominence=%s_width=%s_distance=%s'%(len(peaks),prominence,width,distance)+'.jpg')
    plt.close()  
    return peaks,list_peaks,poperties

def findpeaks2(list_erg,list_eff,peaks,width=50,peaklocateerror=0.005):
    peaks_new=[];porperties_new=[]
    for peak in peaks:
        array_eff=np.array(list_eff)[max(peak-width,0):min(peak+width+1,len(list(list_erg)))]
        peaks_i,porperties_i=sgn.find_peaks(array_eff,prominence=100,width=5)
        peaks_i=[peak_i+max(peak-width,0) for peak_i in peaks_i]
        peaks_new.extend(list(peaks_i))
        # porperties_new=porperties_new+porperties_i
    list_peaks_new=[]
    width2=100
    for peak in peaks_new:
        list_locatedpeaks=[]
        xleft=max(peak-width2,0);xright=min(peak+width2+1,len(list(list_erg)))
        plt.plot(list_erg[xleft:xright],list_eff[xleft:xright])
        plt.axvline(list_erg[peak],linestyle='--')
        for libpeak in lib:
            error=abs((list_erg[peak]/libpeak['Energy']-1))
            if error<peaklocateerror:
                list_locatedpeaks.append(libpeak) 
        list_peaks_new.append(list_locatedpeaks)
        plt.savefig(foldpath+'\\find_peaks\\peak=%s_energy=%.2f.jpg'%(peak,list_erg[peak]))
        plt.close()  
    return peaks_new,porperties_new

# 计算峰面积
# 半高宽拟合
# 使用keV单位
def fwhmfit(erg):
    a=-0.00139;b=0.00517;c=-0.376;erg=erg/1000        
    return (a+b*pow(erg+c*erg**2,0.5))*1000
# 高斯函数拟合
def gauss(x,x0,sigma,h):
    y=h*np.exp(-(x-x0)**2/(2*sigma**2))
    return y
def base(x,a,b):
    y=a*x+b
    return y
def spectrumcompose(x,x0,sigma,h,a,b):
    y=h*np.exp(-(x-x0)**2/(2*sigma**2))+a*x+b
    return y
# 本底+高斯峰拟合
def spectrumfit(list_erg,list_eff,peak,width=50):
    xleft=max(peak-width,0)
    xright=min(peak+width,len(list_eff))
    yleft=list_eff[xleft]
    yright=list_eff[xright]
    a0=(yright-yleft)/(xright-xleft);b0=yleft-xleft*a0
    erg0=list_erg[peak]
    sigma0=abs(fwhmfit(erg0)/2.35)
    amplitude0=list_eff[peak]-b0 
    list_erg=list_erg[xleft:xright]
    list_eff=list_eff[xleft:xright]
    # print(print(erg0,sigma0,amplitude0,a0,b0))
    # print(erg,sigma,amplitude,a,b)
    try:
        erg,sigma,amplitude,a,b=opt.curve_fit(spectrumcompose,list_erg,list_eff,p0=[erg0,sigma0,amplitude0,a0,b0])[0]
    except:
        erg,sigma,amplitude,a,b=erg0,sigma0,amplitude0,a0,b0
    list_fit=[gauss(x,erg,sigma,amplitude)+base(x,a,b) for x in list_erg]
    plt.plot(list_erg,list_eff)
    plt.savefig(foldpath+'\\spectrum_fit\\Unfitted_peak=%s.jpg'%(peak))
    plt.plot(list_erg,list_fit)
    plt.savefig(foldpath+'\\spectrum_fit\\fitted_peak=%s.jpg'%(peak))
    plt.close()
    return erg,sigma,amplitude,a,b
# 峰面积计算方法
# 常数本底扣除+总峰面积法
def areaSum(list_eff,peak,N1=20,N2=5,N3=10):
    counts=0;sigma=0
    for n in range(peak-N3,peak+N3+1):
        counts+=list_eff[n]
        sigma+=list_eff[n]
    background=(sum(list_eff[peak+N1:peak+N1+N2])/N2+sum(list_eff[peak-N1-N2:peak-N1])/N2)/2
    counts-=background*(2*N3+1)
    sigma+=(2*N3+1)*background
    return counts,sigma
# 使用拟合得到的高斯函数系数计算面积
def areaGauss(amplitude,sigma):
    counts=amplitude*pow(2*np.pi*sigma,0.5)/0.29317
    error=sigma
    return counts,error
# def areaGauss(list_erg,list_eff,peak):
    # erg,sigma,amplitude,a,b=spectrumfit(list_erg,list_eff,peak,width=50)
# area2.append(areaGauss(list_erg,list_eff,peak))
# 总峰面积法
def areaTPA(list_eff,peak):
    try:
        list_left=list_eff[peak-50:peak+1]
        list_right=list_eff[peak:peak+51]
        min_left=sgn.argrelmin(np.array(list_left))
        min_right=sgn.argrelmin(np.array(list_right))
        M=min(len(min_left[0]),len(min_right[0]))
        counts=0;sigma=0
        for m in range(M):
            peakleft=min_left[0][-m-1]+peak-50
            peakright=min_right[0][m]+peak
            for n in range(peakleft,peakright+1):
                counts+=list_eff[n]
                sigma+=list_eff[n]
            background=(list_eff[peakleft]+list_eff[peakright])*(peakright-peakleft+1)/2
            counts-=background
            sigma+=(peakright-peakleft+1)*(peakright-peakleft)*(list_eff[peakleft]+list_eff[peakright])/4
        counts=counts/M
        sigma=sigma/M
    except:
        counts,sigma=0,0
    return counts,sigma
# 平均Cowell
def areaCowell(list_erg,list_eff,peak,sigma):
    # erg,sigma,amplitude,a,b=spectrumfit(list_erg,list_eff,peak,width=50)
    Nmin=int(sigma*2.35/0.29317*1)
    Nmax=int(sigma*2.35/0.29317*3)
    counts=0;sigma=0
    for N in range(Nmin,Nmax):
        for n in range(peak-N,peak+N+1):
            counts+=list_eff[n]
            if (n!=(peak-N))and(n!=(peak+N)):
                sigma+=list_eff[n]
        background=(list_eff[peak-N]+list_eff[peak+N])*(2*N+1)/2
        counts-=background
        sigma+=(2*N-1)**2/4*(list_eff[peak-N]+list_eff[peak+N])
        # a=(list_eff[peak-N]-list_eff[peak+N])/(-2*N)
        # b=list_eff[peak-N]-a*(peak-N)
        # list_back=[a*n+b for n in range(len(list_erg))]
        # plt.axvline(x=list_erg[peak],linestyle='--')
        # plt.axvline(x=list_erg[peak-N],linestyle='--')
        # plt.axvline(x=list_erg[peak+N],linestyle='--')
        # plt.plot(list_erg[peak-50:peak+50],list_eff[peak-50:peak+50])
        # plt.plot(list_erg[peak-50:peak+50],list_back[peak-50:peak+50])
        # plt.show()
    counts=counts/(Nmax-Nmin)
    sigma=sigma/(Nmax-Nmin)
    return counts,sigma
# 平均瓦森法找峰
def areaWasson(lis_erg,list_eff,peak,sigma):
    list_left=list_eff[peak-50:peak+1]
    list_right=list_eff[peak:peak+51]
    min_left=sgn.argrelmin(np.array(list_left))
    min_right=sgn.argrelmin(np.array(list_right))
    M=min(len(min_left[0]),len(min_right[0]))
    Nmin=int(sigma*2.35/0.29317)+1
    Nmax=int(sigma*2.35/0.29317*3)
    counts=0;sigma=0;ind=0
    for m in range(M):
        Nmin1=min(Nmin,m);Nmax1=min(Nmax,m)
        for N in range(Nmin1,Nmax1):
            for n in range(peak-N,peak+N+1):
                counts+=list_eff[n]
                if (n!=(peak-N))and(n!=(peak+N)):
                    sigma+=list_eff[n]
            background=(list_eff[peak-N]+list_eff[peak+N])*(2*N+1)/2
            counts-=background
            sigma+=(2*N-1)**2/4*(list_eff[peak-N]+list_eff[peak+N])
            ind+=1
    counts=counts/ind
    sigma=sigma/ind
    return counts,sigma

areaMethodlist=['Sum','Gauss','TPA','Cowell','Wasson']
def areaMethod(list_erg,list_eff,peak,amplitude,sigma,method='Sum'):
    if method=='Sum':
        area,error=areaSum(list_eff,peak)
    elif method=='Gauss':
        area,error=areaGauss(amplitude,peak)
    elif method=='TPA':
        area,error=areaTPA(list_eff,peak)
    elif method=='Cowell':
        area,error=areaCowell(list_erg,list_eff,peak,sigma)
    elif method=='Wassson':
        area,error=areaWasson(list_erg,list_eff,peak,sigma)
    else:
        area,error=0,0
    return area,error
# 使用findpeaks
os.system(r'cd C:\Users\Alber\AppData\Local\VirtualStore\Program Files (x86)\LANL\MCNP5\全能谱分析 & echo y|del .\find_peaks\*.*')
os.system(r'cd C:\Users\Alber\AppData\Local\VirtualStore\Program Files (x86)\LANL\MCNP5\全能谱分析 & echo y|del .\spectrum_fit\*.*')
list_erg,list_eff=spectrumread(filepath)
lib=loadlib(libpath)
peaks,list_peaks,poperties=findpeaks(list_erg,list_eff)
peaks,poperties=findpeaks2(list_erg,list_eff,peaks)
x=1
for peak in peaks:
    spectrumfit(list_erg,list_eff,peak)
    

# for peak in peaks:
#     # peak=peaks[4]
#     erg,sigma,amplitude,a,b=spectrumfit(list_erg,list_eff,peak)
#     list_fit=[gauss(x,erg,sigma,amplitude)+base(x,a,b) for x in list_erg]
#     # print('peak=',peak,'erg=',erg,'sigma=',sigma,'amplitude=',amplitude,a,b,'\n')
#     list_area=[];list_error=[]
#     for method in areaMethodlist:
#         area,error=areaMethod(list_erg,list_eff,peak,amplitude,sigma,method)
#         list_area.append(area)
#         list_error.append(error)
#     print(list_area)
# plt.plot(peaks,list_area,label=method)
# plt.errorbar(peaks,list_area,list_error)
# plt.legend()
# plt.show()

# area1,error1=areaSum(list_eff,peak)
# area2,error2=areaTPA(list_eff,peak)
# area3,error3=areaGauss(list_eff,peak)
# area4,error4=areaCowell(list_erg,list_eff,peak,sigma)
# area5,error5=areaWasson(list_erg,list_eff,peak,sigma)


# plt.plot(list_erg[peak-50:peak+50],list_eff[peak-50:peak+50])
# plt.plot(list_erg[peak-50:peak+50],list_fit[peak-50:peak+50])
# plt.savefig(foldpath+str(peak)+'.jpg')
# plt.close()
# area1=[];area2=[];area3=[];area4=[]
# for i in range(len(peaks)):
#     peak=peaks[i]
#     erg,sigma,amplitude,a,b=spectrumfit(list_erg,list_eff,peak)
#     list_fit=[gauss(x,erg,sigma,amplitude)+base(x,a,b) for x in list_erg]
#     print(peak,erg,sigma,amplitude,a,b,'\n')
#     plt.plot(list_erg[peak-50:peak+50],list_eff[peak-50:peak+50])
#     plt.plot(list_erg[peak-50:peak+50],list_fit[peak-50:peak+50])
#     plt.savefig(foldpath+str(peak)+'kev.jpg')
#     plt.close()
    # area1.append(areaSum(list_eff,peak,N1=20,N2=5,N3=10))
    # area2.append(areaGauss(amplitude,sigma))
    # # area2.append(areaGauss(list_erg,list_eff,peak))
    # area3.append(areaTPA(list_eff,peak)[0])
    # area4.append(areaCowell(list_erg,list_eff,peak,sigma)[0])
    # print(area1,area2,area3,area4)


# plt.plot(peaks,area1,label='Counts')
# plt.plot(peaks,area2,label='Gauss')
# plt.plot(peaks,area3,label='TPA')
# plt.plot(peaks,area4,label='Cowell')
# plt.legend()
# plt.ylim(1e-6,1e-3)
# plt.show()

# axe=plt.subplot(111)
# spectrumplot(axe,list_erg,list_eff,xlogscale=1,ylogscale=2)
# spectrumplot(axe,np.array(list_erg)[peaks],np.array(list_eff)[peaks],para='ro',xlogscale=1,ylogscale=2)
# axe.set_ylim(1e-7,1e-2)
# plt.show()