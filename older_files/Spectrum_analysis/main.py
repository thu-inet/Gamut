import os
import matplotlib.pyplot as plt
import peaksearching as ps
import nuclib as lib
import spectrum as sp

Cs137 = sp.Spectrum()
Cs137.import_GV('E:\\Spectrum_analysis\\Spectrum_data\\GM20燃耗球\\2020-07-21-1800s-026-目标球（反）.chn')
# Cs137.ergcal.recalibrate(list_chn=[2247,5058],list_erg=[661.66,1489.76],L=8192)
# Cs137.import_GV(r'C:\Users\alber\Desktop\WJ.Spe')
# Cs137.smooth_fourier(factor=0.3)
# Cs137.smooth_wavelet(wavelet='sym8',level=3)

Cs137.search_simple()
Cs137.add_peak(386,393)
Cs137.add_peak(1444,1462)
# Cs137.add_peak(1570,1580)
Cs137.add_peak(4091,4100)
Cs137.singlematch(1206,NYlibpath='E:\\Spectrum_Analysis\\Spectrum_data\\NYlib.xml.out',tolerance=0.005)
# print(Cs137.areaTPA(erg=1204,avg_number=2,width_base=2,base_method=0))
# Cs137.draw_spectrum(NYlibpath='E:\\Spectrum_Analysis\\Spectrum_data\\NYlib.xml.out')
# Cs137.report(NYlibpath='E:\\Spectrum_Analysis\\Spectrum_data\\NYlib.xml.out')



# 比较小波光滑效果
# Cs137.smooth_fourier(factor=0.3)
# fig=plt.figure()
# for i in range(1,4):
#     ax=plt.subplot(310+i)
#     # Cs137.smooth_wavelet(wavelet='sym3',level=i)
#     Cs137.smooth_fourier(factor=0.2*i)
    
#     ax.plot(Cs137.list_eff_unsmoothed[1800:2100],'b.',label='O')
#     ax.plot(Cs137.list_eff[1800:2100],label=str(i))
#     plt.yscale('log')
#     plt.legend()
# plt.show()

# plt.plot(Cs137.list_erg[2400:2800],Cs137.list_eff_unsmoothed[2400:2800])
# plt.plot(Cs137.list_erg[2400:2800],Cs137.list_eff[2400:2800])
# plt.plot(Cs137.list_erg[2000:2800],Cs137.list_eff_unsmoothed[2000:2800])
# plt.plot(Cs137.list_erg[2000:2800],Cs137.list_background[2000:2800])
# plt.yscale('log')
# plt.show()


# Cs137.SNIP()
# plt.plot(Cs137.list_erg,Cs137.list_eff_unsmoothed)
# plt.plot(Cs137.list_erg,Cs137.list_background)
# plt.yscale('log')
# plt.show()

# GVlib=nuclib.GVlib_set('D:\\User\\Reports\\Suspect.txt')
# # effcalx=effcal.Effcal_set(c0=-5.391662E-001,c1=2.947627E-001,c2=-4.814516E-008)
# effcalx=effcal.Effcal_set(c0=0.000000E+000,c1=2.931768E-001,c2=0.000000E+000)
# filepath=os.path.abspath('.')
# # list_erg,list_eff=read.read(filepath+'\\Spectrum_data\\GM20燃耗球\\2020-07-21-1800s-023-目标球（反）.Chn')
# list_erg,list_eff=read.read(r'E:\\Spectrum_analysis\\Spectrum_data\\GEM20效率刻度能谱\\Eu-152.spe')
# # list_erg,list_eff=read.read('C:\\Users\\Alber\\Desktop\\WJ.spe')
# list_eff=smoothing.savgol2(list_eff,width=9,order=3)
# peaks=ps.simpleseach(list_eff,prominence=1.1,halfwidth=3)
# file=open(filepath+'\\rpt.txt','w')
# file.write('Nuclide  Nuc_Energy(keV) Peak_Energy(keV) Left_chn Right_chn \t Total Area \t Uncertainty \n')
# # for peak in peaks:
# #     peak.energy=effcalx.chn2erg(peak.centroid)
# # for peak in peaks[:-2]:
# #     area,uncertainty=areaTPA.areaTPA(list_erg,list_eff,peak=peak.centroid)[0:2]
# #     nuc=GVlib.peak_match([peak])[0]
# #     if area>100:
# #         if nuc == 'Unknown':
# #             file.write('Unknown \t          \t %-8.2f \t\t %-8i \t %-8i \t %-10.2f \t %-8.4f\n'%(peak.energy,peak.left,peak.right,area,uncertainty))
# #         else:
# #             file.write('%-8s \t %-8s \t %-8.2f \t\t %-8i \t %-8i \t %-10.2f \t %-8.4f\n'%(nuc.nuclide,nuc.energy,peak.energy,peak.left,peak.right,area,uncertainty))
# # file.close()

# print(areaTPA.areaTPA(list_erg,list_eff,peak=3701))
# Cs137=Spectrum()