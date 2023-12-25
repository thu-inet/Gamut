import spectrumread
import numpy as np
import matplotlib.pyplot as plt
import numpy.fft 
import scipy.signal as sgn 
import xlwt
import scipy.optimize as op
pi=3.1415926

filename='Eu-152.spe'
filepath=r'C:\Users\Alber\AppData\Local\VirtualStore\Program Files (x86)\LANL\MCNP5\全能谱分析\\'
filepath=filepath+filename
list_erg,list_eff=spectrumread.spectrumread(filepath)

# smooth by barycenter method
# def smoothBarycenter(list):
#  smooth by fourier transformation
def smoothSavgol(list_eff,WindowWidth=13,PolynomeOrder=5):
    list_eff_smoothed=sgn.savgol_filter(list_eff,WindowWidth,PolynomeOrder)
    return list_eff_smoothed

#  smooth by fourier transformation
def weight(i,sigma=0.6):
    A=1
    omiga=(4092-abs(i-4092))/8192*2*pi
    y=A*np.exp(-0.5*sigma**2*omiga**2)
    return y

# def smoothFourier(list_eff,weightfunction=weight):
#     list_eff_transferred=numpy.fft.fft(list_eff)
#     list_weight=[weightfunction(i) for i in range(len(list_eff))]
#     list_eff_transferred_weighted=[list_eff_transferred[i]*list_weight[i] for i in range(len(list_eff))]
#     list_eff_weighted_retransferred=numpy.fft.ifft(list_eff_transferred_weighted)
#     return list_eff_weighted_retransferred

list_eff2=smoothSavgol(list_eff)
list_eff_transferred=numpy.fft.fft(list_eff)
list_weight=[weight(i) for i in range(len(list_eff))]
list_eff_transferred_weighted=[list_eff_transferred[i]*list_weight[i] for i in range(len(list_eff))]
list_eff_weighted_retransferred=numpy.fft.ifft(list_eff_transferred_weighted)

# list_eff3=smoothFourier(list_eff)
# plt.plot(list_eff[600:1000],'b-')
# plt.plot(list_eff2[600:1000],'r-')
# plt.plot(list_eff3[600:1000],'k-')
# plt.yscale('log')

plt.plot(list_eff[1000:2000],'b')
# plt.plot(list_weight)
# plt.plot(list_eff_transferred,'r.')
# plt.plot(list_eff_transferred_weighted,'k')
plt.plot(list_eff_weighted_retransferred[1000:2000],'c')
plt.show()

# fig=plt.figure()
# axe1=plt.subplot(111)

# axe1.plot(erg,cts,linewidth=0.5,label='origin')
# axe1.plot(erg,cts_f,'b',label='Spc',linewidth=0.5)
# axe1.plot(np.linspace(0,2*pi,1500),cts_f,'b',label='Spc',linewidth=0.5)


# book=xlwt.Workbook(encoding='utf-8')
# for sigma in [0.1,0.2,0.4,0.6,0.8,1,1.2,1.4,1.6,1.8,2,4]:
#     cts_filter=[filter(i,sigma) for i in erg]
#     cts_f_filter=[cts_f[i]*cts_filter[i] for i in range(len(cts_f))]
#     cts_new=f.ifft(cts_f_filter)
#     plt.plot(erg,cts,linewidth=0.5,label='origin')
#     plt.plot(erg,cts_new,linewidth=1,label='sigma='+str(sigma)[:6])
#     h=2
#     while True:
#         peak_index,_=sgn.find_peaks(cts_new,height=h)
#         if len(peak_index)>50:
#             h=h*2
#         elif len(peak_index)<10:
#             h=h/10
#         else:
#             break
#     sheet=book.add_sheet('sigma='+str(sigma)[:6])
#     sheet.write(0,0,'channel')
#     sheet.write(0,1,'energy')
#     sheet.write(0,2,'efficiency')
#     sheet.write(0,3,'FWHM')
#     for i in range(1500):
#         list_erg[i]=list_erg[i]/1500*1.5
#     list_peak_erg=[]
#     list_peak_cts=[]
#     list_peak_fwhm=[]
#     k=0 
#     for peak in peak_index:
#         m=max(peak-30,0)
#         M=min(peak+30,len(list_erg))
#         ergi=list_erg[m:M]
#         ctsi=cts_new[m:M]
#         x0,sigma,h,base=op.curve_fit(Gauss,ergi,ctsi,p0=[list_erg[peak],2e-3,list_cts[peak],0])[0]
#         print(peak,'erg=',x0,'sigma=',sigma,'eff=',h,'\t')
#         if x0>0.2:
#             k+=1
#             try:
#                 for x in range(len(list_peak)):
#                     if (x0-list_peak[x])/x0<0.05:
#                         eff_i=list_i[x]
#                 eff=h*pow(2*pi,0.5)*sigma*sum_i/eff_i*0.015/2/0.29317e-3
#                 fwhm=2.335*sigma
#                 list_peak_erg.append(x0)
#                 list_peak_cts.append(eff)
#                 list_peak_fwhm.append(sigma*2.335)
#                 sheet.write(k,0,str(peak))
#                 sheet.write(k,1,str(x0))
#                 sheet.write(k,2,str(eff))
#                 sheet.write(k,3,str(fwhm))
#                 sheet.write(k,4,'erg='+str(x0)+'sigma='+str(sigma)+'h='+str(h)+'\t')
#             except:
#                 x=1
#     book.save(filepath+'smoothing detection efficiency'+'.xls')


#     plt.legend()
#     plt.savefig(filepath+'sigma='+str(sigma)[:6]+'.jpg')
#     plt.close()



# plt.ylim(1.5,2.5)



    
    
  