import spectrum
import os
import matplotlib.pyplot as plt
import xlwt

# parameter settings
spectrum_filename='Eu-152'
time_ratio=230335/62139
time_ratio_square=time_ratio**2
filepath=r'C:\Users\Alber\AppData\Local\VirtualStore\Program Files (x86)\LANL\MCNP5\全能谱分析\Data'
list_erg,list_area_baseline=spectrum.read(filepath+'\\本底.spe')
list_erg,list_area_total=spectrum.read(filepath+'\\%s.spe'%(spectrum_filename))
list_peak_channel=[416,837,1176,2656,3288,3704,3792,4803]
peak_channel=3288
report=open('AreaUncertaintyReduction.txt','w',encoding='utf-8')

list_area_net=[]
for channel in range(len(list_erg)):
    list_area_net.append(list_area_total[channel]-time_ratio*list_area_baseline[channel])
spectrum.write(list_area_net,filepath='Eu-152-no-baseline.spe')


# 使用两个能谱做差
report.write('# 使用两个能谱做差\n')
report.write(' roi_width_r=    -10-    -11-    -12-    -13-    -14-    -15-    -16-    -17-     -18-    -19-    -20-\n')
for roi_width_l in range(10,21):
    report.write('roi_width_l=%i\t'%(roi_width_l))
    for roi_width_r in range(10,21):
        area=0
        variance=0
        for channel in range(peak_channel-roi_width_l,peak_channel+roi_width_r+1):
            area+=list_area_total[channel]-time_ratio*list_area_baseline[channel]
            variance+=list_area_total[channel]+time_ratio_square*list_area_baseline[channel]
        std_variance=variance**0.5
        uncertainty=std_variance/area
        # report.write('Width_left=%i  Width_right=%i  Area=%6f  Error=%6f   Sigma=%3.3f  \n'%(roi_width_l,roi_width_r,area,std_variance,uncertainty*100))
        report.write('%3.3f%% \t'%(uncertainty*100))
    report.write('\n')
report.write('\n')

report.write('Width_left    Width_right   Net Area     Error    Sigma\n')
for roi_width_l in range(10,21):
    for roi_width_r in range(10,21):
        area=0
        variance=0
        for channel in range(peak_channel-roi_width_l,peak_channel+roi_width_r+1):
            area+=list_area_total[channel]-time_ratio*list_area_baseline[channel]
            variance+=list_area_total[channel]+time_ratio_square*list_area_baseline[channel]
        std_variance=variance**0.5
        uncertainty=std_variance/area
        if roi_width_r==10:
            report.write('%-10.4s    %-10.4s    %-8.8s    %-6.6s   %-6.6s  \n'%(roi_width_l,roi_width_r,area,std_variance,uncertainty*100))
        else:
            report.write('              %-10.4s    %-8.8s    %-6.6s   %-6.6s  \n'%(roi_width_r,area,std_variance,uncertainty*100))


# 使用一个能谱扣除本底
roi_width_r=roi_width_l=40
report.write('# 使用一个能谱扣除本底\n')
for roi_base in range(5,30):
    base_left=sum(list_area_total[peak_channel-roi_width_l-roi_base:peak_channel-roi_width_r])/roi_base
    base_right=sum(list_area_total[peak_channel+roi_width_l-roi_base:peak_channel+roi_width_r])/roi_base
    area=0;variance=0
    for channel in range(peak_channel-roi_width_l,peak_channel+roi_width_r+1):
        area+=list_area_total[channel]
        variance+=list_area_total[channel]
    area-=(base_left+base_right)/2*(roi_width_l+roi_width_r+1)
    variance+=(base_left+base_right)/roi_base/4*(roi_width_l+roi_width_r+1)**2
    std_variance=variance**0.5
    uncertainty=std_variance/area
    report.write('ROi_base=%4i   Area=%4s  Error=%4s   Sigma=%4s  \n'%(roi_base,area,std_variance,uncertainty*100))


report.close()

# list_erg,total_counts=read(r'C:\Users\Alber\AppData\Local\VirtualStore\Program Files (x86)\LANL\MCNP5\全能谱分析\Data'+'\\Cs-137.spe')
# list_erg,total_counts=read(r'C:\Users\Alber\AppData\Local\VirtualStore\Program Files (x86)\LANL\MCNP5\全能谱分析\Data'+'\\Co-60.spe')
list_erg,total_counts=read(r'C:\Users\Alber\AppData\Local\VirtualStore\Program Files (x86)\LANL\MCNP5\全能谱分析\Data'+'\\Eu-152.spe')

list_peak_channel=[416,837,1176,2656,3288,3704,3792,4803]
# print(str(areaTPA(total_counts,peak=3288,avg_number=2,width_base=2)))

report=open('report.txt','w')
for peaks in list_peak_channel:
    print(peaks)
    report.write('peak=%i\t'%(peaks)+str(areaTPA(total_counts,peak=peaks,avg_number=3,width_base=6))+'\n')
report.close()
