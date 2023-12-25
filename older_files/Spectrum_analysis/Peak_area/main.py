import spectrum
import matplotlib.pyplot as plt
import numpy as np
import scipy.optimize as opt
import sko
import time
# list_erg,total_counts=read(r'C:\Users\Alber\AppData\Local\VirtualStore\Program Files (x86)\LANL\MCNP5\全能谱分析\Data'+'\\Cs-137.spe')
# list_erg,total_counts=read(r'C:\Users\Alber\AppData\Local\VirtualStore\Program Files (x86)\LANL\MCNP5\全能谱分析\Data'+'\\Co-60.spe')
# list_erg,total_counts=read(r'C:\Users\Alber\AppData\Local\VirtualStore\Program Files (x86)\LANL\MCNP5\全能谱分析\Data'+'\\Eu-152.spe')
# list_erg,total_counts=spectrum.read(r'E:\\MCNP_data\\spectrum_analysis\Data'+'\\2020-07-21-1800s-023-目标球.spe')
# list_erg,base_counts=spectrum.read(r'E:\\MCNP_data\\spectrum_analysis\Data'+'\\Base.spe')

# peak=2053
# le=list_erg[peak-30:peak+30]
# lf=total_counts[peak-30:peak+30]
# list_erg=le
# list_eff=lf
# global id
# id=0
# def triple_gausspeak(p):
#     e0,h0,s0,e1,h1,s1,e2,h2,s2,b,d=p
#     sum_error=0
#     global id
#     list_y=[]
#     for i in range(len(list_erg)):
#         erg=list_erg[i]
#         eff=list_eff[i]
#         y=h0*np.exp(-(erg-e0)**2/2/s0**2)+h1*np.exp(-(erg-e1)**2/2/s1**2)+h2*np.exp(-(erg-e2)**2/2/s2**2)+b*np.exp(-erg/d)
#         list_y.append(y)
#         sum_error+=(y/eff-1)**2
#     plt.plot(le,lf)
#     plt.plot(le,list_y)
#     plt.savefig("E:\\mcnp_data\\x\\%8s.jpg"%(id))
#     plt.ylim([0,4000])
#     plt.close()
#     id+=1
#     # print(sum_error)
#     return sum_error
# def gausspeak(p):
#     e0,h0,s0,b,d=p
#     sum_error=0
#     global id
#     list_y=[]
#     for i in range(len(list_erg)):
#         erg=list_erg[i]
#         eff=list_eff[i]
#         # y=h0*np.exp(-(erg-e0)**2/s0)+b*(erg-600)+d
#         y=h0*np.exp(-(erg-e0)**2/s0)+b*(erg-600)+d
#         list_y.append(y)
#         sum_error+=(y/eff-1)**2
#     plt.plot(le,lf)
#     plt.plot(le,list_y)
#     plt.savefig("E:\\mcnp_data\\x\\%8s.jpg"%(id))
#     plt.ylim([0,4000])
#     plt.close()
#     id+=1
#     # print(sum_error)
#     return sum_error
# GArun=sko.GA.GA(gausspeak,n_dim=5,size_pop=200,max_iter=400,\
#     prob_mut=0.001,lb=[600,1e2,1,-10,1800],ub=[605,1e4,100,10,2200],precision=1e2)
# # GEr,GEh,CUr,CUh,SDLd,TDLd,Ald,XD,YD=P
# bestx,besty=GArun.run()
# print('------------------------------------')
# print(bestx,besty)



# param=opt.curve_fit(triple_gausspeak,le,lf,p0=[605,1,1000,1,1,601,1,1,610,1,1])[0]
# print(param)
# while True:
#     param=list(param)
#     param=opt.curve_fit(triple_gausspeak,le,lf,p0=param)[0]
#     print(param)
# peak=2053
# le=list_erg[peak-100:peak-20]
# le.extend(list_erg[peak+20:peak+100])
# lf=total_counts[peak-100:peak-20]
# lf.extend(total_counts[peak+20:peak+100])
# plt.plot(le,lf)

# plt.show()
# peak=2702
# w=200
# list_eff_back=spectrum.SNIP(list_erg,total_counts)
# plt.plot(list_erg[peak-w:peak+w],total_counts[peak-w:peak+w])
# plt.plot(list_erg[peak-w:peak+w],list_eff_back[peak-w:peak+w])
# plt.show()
# print(str(spectrum.areaTPA(list_erg,total_counts,peak=2600,avg_number=5,width_base=5,base_method=1)))  
# for m in range(1,30):
#     for n in range(1,30):
        # print(m,n,str(spectrum.areaTPA(list_erg,total_counts,peak=2246,avg_number=m,width_base=n)))  
x=1
# list_peak_channel=[416,837,1176,2656,3288,3704,3792,4803]
# for peaks in list_peak_channel:
#     print(peaks)
#     report.write('peak=%i\t'%(peaks)+str(areaTPA(total_counts,peak=peaks,avg_number=3,width_base=6))+'\n')
# report.close()
# lst=[]
# for i in range(len(total_counts)):
#     lst.append(total_counts[i]-base_counts[i]*230335/62139)
# plt.plot(list_erg[3650:3750],lst[3650:3750])
# plt.show()
# for l in range(2645,2649):
#     for r in range(2664,2670):
#         s=sum(total_counts[l:r+1])-sum(base_counts[l:r+1])*230335/62139
#         print(l,r,sum(total_counts[l:r+1]),sum(base_counts[l:r+1]),s)    

# list_peak=[112.42,134.45,428.31,497.25,511.98,604.58,621.87,661.37,696.12,723.74,756.19,765.24,795.2,872.67,1049.11,1203.05,1458.89,1486.87,1631.05,2181.75]

# def find(list_erg,peak):
#     l=len(list_erg)
#     i1=0
#     i2=l     
#     while 1:
#         i=int((i1+i2)/2)
#         erg=list_erg[i]
#         if abs(erg/peak-1)<=1.5e-3:
#             return i
#         elif erg<peak:
#             i1=i
#         else:
#             i2=i

# for peak in list_peak:
#     ind=find(list_erg,peak)
#     print(peak,ind,str(areaTPA(total_counts,peak=ind,avg_number=5,width_base=10)))




