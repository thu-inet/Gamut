#输入一个mcnp的输出文件，给出对应的能谱
import numpy as np
import matplotlib.pyplot as plt
import numpy.fft as f
import scipy.signal as sgn 
import xlwt
import scipy.optimize as op
pi=3.1415926
filename='o0_0-051905'
filepath=r'C:\Users\Alber\AppData\Local\VirtualStore\Program Files (x86)\LANL\MCNP5\全能谱分析\\'
filepath=filepath+filename
fileopen=open(filepath)
fileread=fileopen.read()
list_erg=[]
list_cts=[]
list_peak=[0.1217817,0.24469758,0.344278573,0.661657,0.778904,0.964079,1.0858682474,1.1120744,1.173228,1.332492,1.4080063]
list_i=[0.28256,0.07583,0.2654,0.851,0.12942,0.14605,0.10207,0.13644,0.9985,0.999826,0.21005]
sum_i=sum(list_i)
def Gauss(x,x0,sigma,h,base):
    return h*np.exp(-(x-x0)**2/sigma**2/2)+base
# def noise():
#     n=abs(np.random.normal(0,1))
#     return n

ind_pre=fileread.index('1tally   8        nps =')
ind=fileread.index('0.0000E+00',ind_pre)
sum_eff=0
erg=0
while True:
    try:
        erg+=1
        cts=int(float(fileread[ind+13:ind+24])*10000000)
        cts=np.log10(cts+1)
        list_erg.append(erg)
        list_cts.append(cts)
        li=len(str(cts))
        ind+=36
    except:
        del list_erg[-1]
        del list_cts[-1]
        break

for i in range(10):
    list_cts[i]=list_cts[10]


erg=list_erg;cts=list_cts
# erg=list_erg[500:1000]
# cts=list_cts[500:1000]



# plt.plot(erg,cts,'g',label='Spctrum',linewidth=0.5)


# # smooth by solvog
# cts_s=sgn.savgol_filter(cts,5,2)
# plt.plot(erg,cts_s,'b',label='cts_s_52',linewidth=0.5)
# cts_s=sgn.savgol_filter(cts,7,3)
# plt.plot(erg,cts_s,'y',label='cts_s_73',linewidth=0.5)
# cts_s=sgn.savgol_filter(cts,25,8)
# plt.plot(erg,cts_s,'r',label='cts_s_94',linewidth=0.5)

# smooth by fourier transformation

cts_f=f.fft(cts)

# def filter(i):
#     A=1
#     a=-0.00139;b=0.00517;c=-0.376
#     erg=i*0.29317*0.001
#     H=a+b*pow(erg+erg**2*c,0.5)
#     H=H*1000/0.29317
#     sigma=H/2.355
#     omiga=(i-500)/500*2*pi
#     f=A*np.exp(-0.5*sigma**2*omiga**2)
#     return f
# cts_filter=[filter(i) for i in erg]
# cts_f_filter=[cts_f[i]*cts_filter[i] for i in range(len(cts_f))]
# cts_new=f.ifft(cts_f_filter)

# fig=plt.figure()
# axe1=plt.subplot(221)

# axe1.plot(erg,cts,'b',label='Spc',linewidth=0.5)
# axe1.plot(erg,cts_new,'r',linewidth=0.5,label='Spc_filter')

# axe2=plt.subplot(222)
# axe3=plt.subplot(223)
# axe4=plt.subplot(224)
# axe2.plot(erg,cts_f,'r',label='能频谱')
# axe3.plot(erg,cts_filter,linewidth=0.5)
# axe4.plot(erg,cts_f_filter,linewidth=0.5)


def filter(i,sigma):
    A=1
    # omiga=(i-500)/500*2*pi
    # omiga=i/1500*2*pi
    omiga2=(750-abs(i-750))/1500*2*pi
    f=A*np.exp(-0.5*sigma**2*omiga2**2)
    return f

# fig=plt.figure()
# axe1=plt.subplot(111)

# axe1.plot(erg,cts,linewidth=0.5,label='origin')
# axe1.plot(erg,cts_f,'b',label='Spc',linewidth=0.5)
# axe1.plot(np.linspace(0,2*pi,1500),cts_f,'b',label='Spc',linewidth=0.5)


book=xlwt.Workbook(encoding='utf-8')
for sigma in [0.1,0.2,0.4,0.6,0.8,1,1.2,1.4,1.6,1.8,2,4]:
    cts_filter=[filter(i,sigma) for i in erg]
    cts_f_filter=[cts_f[i]*cts_filter[i] for i in range(len(cts_f))]
    cts_new=f.ifft(cts_f_filter)
    plt.plot(erg,cts,linewidth=0.5,label='origin')
    plt.plot(erg,cts_new,linewidth=1,label='sigma='+str(sigma)[:6])
    h=2
    while True:
        peak_index,_=sgn.find_peaks(cts_new,height=h)
        if len(peak_index)>50:
            h=h*2
        elif len(peak_index)<10:
            h=h/10
        else:
            break
    sheet=book.add_sheet('sigma='+str(sigma)[:6])
    sheet.write(0,0,'channel')
    sheet.write(0,1,'energy')
    sheet.write(0,2,'efficiency')
    sheet.write(0,3,'FWHM')
    for i in range(1500):
        list_erg[i]=list_erg[i]/1500*1.5
    list_peak_erg=[]
    list_peak_cts=[]
    list_peak_fwhm=[]
    k=0 
    for peak in peak_index:
        m=max(peak-30,0)
        M=min(peak+30,len(list_erg))
        ergi=list_erg[m:M]
        ctsi=cts_new[m:M]
        x0,sigma,h,base=op.curve_fit(Gauss,ergi,ctsi,p0=[list_erg[peak],2e-3,list_cts[peak],0])[0]
        print(peak,'erg=',x0,'sigma=',sigma,'eff=',h,'\t')
        if x0>0.2:
            k+=1
            try:
                for x in range(len(list_peak)):
                    if (x0-list_peak[x])/x0<0.05:
                        eff_i=list_i[x]
                eff=h*pow(2*pi,0.5)*sigma*sum_i/eff_i*0.015/2/0.29317e-3
                fwhm=2.335*sigma
                list_peak_erg.append(x0)
                list_peak_cts.append(eff)
                list_peak_fwhm.append(sigma*2.335)
                sheet.write(k,0,str(peak))
                sheet.write(k,1,str(x0))
                sheet.write(k,2,str(eff))
                sheet.write(k,3,str(fwhm))
                sheet.write(k,4,'erg='+str(x0)+'sigma='+str(sigma)+'h='+str(h)+'\t')
            except:
                x=1
    book.save(filepath+'smoothing detection efficiency'+'.xls')


    plt.legend()
    plt.savefig(filepath+'sigma='+str(sigma)[:6]+'.jpg')
    plt.close()



# plt.ylim(1.5,2.5)



    
    
  