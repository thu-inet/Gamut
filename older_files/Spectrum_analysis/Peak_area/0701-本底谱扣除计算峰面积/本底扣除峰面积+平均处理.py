import spectrum 
import os
import numpy as np
import scipy.optimize as opt
import matplotlib.pyplot as plt
import xlwt

# parameter settings
filename='Eu-152'
peakl=[416,837,1176,2656,3288,3704,3792,4803]
# peakl=[4000,4544]
# peakl=[2256]
for peak in peakl:
    #Co-60 86483,Cs-137:88562,Eu:230335
    r=230335/62139;r2=r**2
    path=os.path.abspath('.')
    erg,total=spectrum.read(path+'\\全能谱分析\%s.spe'%(filename))
    erg,baseline=spectrum.read(path+'\\全能谱分析\\本底.spe')
    net=[total[i]-baseline[i]*r for i in range(len(total))]
    width=100
    
    xlt=xlwt.Workbook(encoding='utf-8')
    baser=[r*base for base in baseline[peak-width:peak+width+1]]
    plt.plot(erg[peak-width:peak+width+1],total[peak-width:peak+width+1],label='Total')
    plt.plot(erg[peak-width:peak+width+1],baser,label='Baseline')
    plt.plot(erg[peak-width:peak+width+1],net[peak-width:peak+width+1],label='Net')
    plt.axhline(y=0,linestyle='--')
    plt.title('Peak=%s Energy=%skeV'%(peak,erg[peak]))
    plt.xlabel('Energy');plt.ylabel('Counts')
    plt.legend()
    plt.savefig(path+'\\全能谱分析\\%s-%s.png'%(filename,peak))
    plt.close()

    # Excel preparation
    sheet=xlt.add_sheet('Spectrum data')
    sheet.write(0,0,'Channel')
    sheet.write(0,1,'Energy')
    sheet.write(0,2,'Total counts')
    sheet.write(0,3,'Baseline counts')
    sheet.write(0,4,'Net counts')
    for i in range(peak-width,peak+width+1):
        j=i-peak+width
        sheet.write(j+1,0,i)
        sheet.write(j+1,1,erg[i])
        sheet.write(j+1,2,total[i])
        sheet.write(j+1,3,baseline[i])
        sheet.write(j+1,4,net[i])

    # method1-peak area counting
    sheet=xlt.add_sheet('Total peak area method')
    sheet.write(0,0,'Half widths/channels')
    sheet.write(0,1,'Area counts')
    sheet.write(0,2,'Sigma')
    sheet.write(0,3,'Uncertainty/%')
    sheet.write(0,4,'Uncertainty(compared to mean area)/%')
    list_area=[]
    list_error=[]
    imin=10;imax=30
    for i0 in range(imin,imax+1):
        area=0;sigma2=0
        for i in range(peak-i0,peak+i0+1):
            area+=net[i]
            sigma2+=total[i]+r2*baseline[i]
        error=pow(sigma2,0.5)/area
        sheet.write(i0-imin+1,0,i0)
        sheet.write(i0-imin+1,1,area)
        sheet.write(i0-imin+1,2,pow(sigma2,0.5))
        sheet.write(i0-imin+1,3,error)
        list_area.append(area)
        list_error.append(error)
    mean=sum(list_area)/(imax-imin+1)
    for i0 in range(imin,imax+1):
        sheet.write(i0-imin+1,4,(list_area[i0-imin]/mean-1)*100)
    # max_sum=max(list_sum)
    # max_error=max(list_error)
    # l=len(list_sum)
    # plt.plot(np.linspace(10,20,l),list_sum)
    # plt.axhline(max_sum*0.995)
    # plt.ylim((30000,40000))
    # plt.plot(np.linspace(10,20,l),list_error)

    # method2-average peak area counting
    sheet=xlt.add_sheet('Avarage total peak area method')
    sheet.write(0,0,'Min half widths/channels')
    sheet.write(0,1,'Max half widths/channels')
    sheet.write(0,2,'Area counts')
    sheet.write(0,3,'Sigma')
    sheet.write(0,4,'Uncertainty/%')
    sheet.write(0,5,'Uncertainty(compared to mean area)/%')
    list_area2=[]
    list_error2=[]
    for i0 in range(imin,imax+1):
        area=sum(list_area[0:i0-imin+1])/(i0-imin+1)
        sigma2=0
        for i in range(peak-i0,peak+i0+1):
            d=i0+1-max(abs(i-peak),imin)
            sigma2+=(total[i]+r2*baseline[i])*d**2
        sigma2=sigma2/(i0-imin+1)**2
        error=pow(sigma2,0.5)/area
        sheet.write(i0-imin+1,0,imin)
        sheet.write(i0-imin+1,1,i0)
        sheet.write(i0-imin+1,2,area)
        sheet.write(i0-imin+1,3,pow(sigma2,0.5))
        sheet.write(i0-imin+1,4,error)
        list_area2.append(area)
        list_error2.append(error)
    mean=sum(list_area2)/(imax-imin+1)
    for i0 in range(imin,imax+1):
        sheet.write(i0-imin+1,5,(list_area2[i0-imin]/mean-1))

    # method3-simgle spectrum
    # sheet=xlt.add_sheet('Avarage total peak area method')
    # sheet.write(0,0,'Min half widths/channels')
    # sheet.write(0,1,'Max half widths/channels')
    # sheet.write(0,2,'Area counts')
    # sheet.write(0,3,'Sigma2')
    # sheet.write(0,4,'Uncertainty/%')
    # sheet.write(0,5,'Uncertainty(compared to mean area)/%')
    # list_area2=[]
    # list_error2=[]
    # for i0 in range(imin,imax+1):
    #     area=sum(list_area[0:i0-imin+1])/(i0-imin+1)
    #     sigma2=0
    #     for i in range(peak-i0,peak+i0+1):
    #         d=i0+1-max(abs(i-peak),imin)
    #         sigma2+=(total[i]+r2*baseline[i])*d**2
    #     sigma2=sigma2/(i0-imin+1)**2
    #     error=pow(sigma2,0.5)/area
    #     sheet.write(i0-imin+1,0,imin)
    #     sheet.write(i0-imin+1,1,i0)
    #     sheet.write(i0-imin+1,2,area)
    #     sheet.write(i0-imin+1,3,pow(sigma2,0.5))
    #     sheet.write(i0-imin+1,4,error)
    #     list_area2.append(area)
    #     list_error2.append(error)
    # mean=sum(list_area2)/(imax-imin+1)
    # for i0 in range(imin,imax+1):
    #     sheet.write(i0-imin+1,5,(list_area2[i0-imin]/mean-1)*100)
    xlt.save(path+'\\全能谱分析\\%s-%s.xls'%(filename,peak))

# # method3 peakfit
# sheet=xlt.add_sheet('Peakfit method')
# sheet.write(0,0,'Channel')
# sheet.write(0,1,'Energy')
# sheet.write(0,2,'Net counts')
# sheet.write(0,3,'Counts fit')
# sheet.write(0,3,'Fir errot/%')
# netfit=net[peak-width:peak+width+1]
# ergfit=erg[peak-width:peak+width+1]


# def f1(x,h1,e1,s1):
#     y=h1*np.exp(-(x-e1)**2/s1**2)
#     return y
# p0=[netfit[width],ergfit[width],0.5]
# para,covar=opt.curve_fit(f1,ergfit,netfit,p0)
# h1,e1,s1=para[0],para[1],para[2]
# fit=[f1(ergfit[i],h1,e1,s1) for i in range(peak-width,peak+width+1)]

# # def f2(x,h1,e1,s1,h2,e2,s2):
# #     y=h1*np.exp(-(x-e1)**2/s1**2)+h2*np.exp(-(x-e2)**2/s2**2)
# #     return y
# # def f2withbase(x,h1,e1,s1,h2,e2,s2,a,b):
# #     y=h1*np.exp(-(x-e1)**2/s1**2)+h2*np.exp(-(x-e2)**2/s2**2)+a*x+b
# #     return y
# # p0=[netfit[width],ergfit[width],0.5]
# # para=opt.curve_fit(f2,ergfit,netfit,p0)
# # h1,e1,s1,h2,e2,s2,a,b=opt.curve_fit(f,erg2,net,p0)[0]
# # para=[h1,e1,s1,h2,e2,s2,a,b]
# plt.plot(ergfit,total[peak-width:peak+width+1],label='Total')
# plt.plot(ergfit,netfit,label='Net')
# plt.legend()
# plt.plot(ergfit,fit,label='Fit')
# plt.savefig(path+'\\全能谱分析\\%s-%s.png'%(filename,peak))

# for i in range(peak-width,peak+width+1):
#     j=i-peak+width
#     sheet.write(j+1,0,i)
#     sheet.write(j+1,1,erg[i])
#     sheet.write(j+1,2,net[i])
#     sheet.write(j+1,3,fit[j])
#     sheet.write(j+1,3,(fit[j]/net[i]-1)*100)

# print(para)
# # plt.plot(erg2,net,label='net')
# # plt.plot(erg2,fit,label='fit')
# # plt.legend()
# # plt.show()
# sum1=sum(net)
# sum2=h1*pow(2*np.pi,0.5)*s1/0.29317


