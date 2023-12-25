import time
import matplotlib.pyplot as plt
import numpy as np



def draw_spectrum(instance_spectrum,emin=0,emax=2400,ymin=1,ymax=10000,NYlibpath=None):
    if not hasattr(instance_spectrum,'list_gamma_match') or (len(instance_spectrum.list_peaks)!=len(instance_spectrum.list_gamma_match)):
        instance_spectrum.match(NYlibpath=NYlibpath)
    plt.figure(figsize=(30,10))
    # plt.plot(instance_spectrum.list_erg,instance_spectrum.list_eff_unsmoothed,'b.',label='Unsmoothed',markersize=0.5)
    plt.plot(instance_spectrum.list_erg,instance_spectrum.list_eff,'b.',markersize=1,label='Smoothed',linewidth=0.5)
    
    for i_peak,peak in enumerate(instance_spectrum.list_peaks):
        centroid=list(instance_spectrum.list_eff[peak.peak_left_boundary:peak.peak_right_boundary+1]).index(max(instance_spectrum.list_eff[peak.peak_left_boundary:peak.peak_right_boundary+1]))+peak.peak_left_boundary
        if instance_spectrum.list_eff[centroid]-instance_spectrum.list_eff[peak.peak_left_boundary]>20:
            plt.text(instance_spectrum.list_erg[centroid],
                     instance_spectrum.list_eff_unsmoothed[centroid],
                     '  <--%6.2fkeV(%s)'%(instance_spectrum.list_erg[peak.centroid],instance_spectrum.list_gamma_match[i_peak].nuclide),
                     rotation=90,ha='center',va='bottom',fontsize=6)
        plt.plot(instance_spectrum.list_erg[peak.peak_left_boundary:peak.peak_right_boundary+1],instance_spectrum.list_eff[peak.peak_left_boundary:peak.peak_right_boundary+1],'r-',linewidth=0.7)
        plt.fill_between(instance_spectrum.list_erg[peak.peak_left_boundary:peak.peak_right_boundary+1],instance_spectrum.list_eff[peak.peak_left_boundary:peak.peak_right_boundary+1],alpha=0.7)
    # plt.axvline(x=list_erg[peak_boundary_left],ymax=list_eff[peak_boundary_left]/ymax,linestyle='--')
    # plt.axvline(x=list_erg[peak_boundary_right],ymax=list_eff[peak_boundary_right]/ymax,linestyle='--')
    plt.xlabel('Energy, keV')
    plt.ylabel('Counts')
    plt.yscale('log')
    plt.title('Peak Region Plot')
    plt.xticks(np.arange(0,2500,50))
    plt.legend()
    plt.xlim([emin,emax])
    ymax=np.exp(int(np.log(max(instance_spectrum.list_eff))+1))
    plt.ylim([ymin,ymax])
    plt.savefig(r'C:\Users\Alber\Desktop\spectrum.svg')
    plt.show()

def strip(instance_spectrum,background,ratio):
    list_eff_strpied=[]
    for total,back in zip(instance_spectrum.list_eff,background.list_eff):
        if total !=0:
            net=max(int(total-back*ratio+1),0)  
        else:
            net=0
        list_eff_strpied.append(net)
    instance_spectrum.list_eff_unsmoothed=list_eff_strpied
    instance_spectrum.list_eff=list_eff_strpied
    return 0


def report(instance_spectrum,GVlibpath='D:\\User\\Reports\\Suspect.txt',NYlibpath=None):
    if (not hasattr(instance_spectrum,'list_gamma_match')) or (NYlibpath!=None):
        instance_spectrum.match(GVlibpath=GVlibpath,NYlibpath=NYlibpath)
    report= 'No.        Channel      Energy/keV         Isotope      Energy/keV    Half_life(yrs)       Intensity      Net_area       Error      Total_area       Peak boundaries  Back boundaries \n'
           # 57           3567         1049.15          Rh-106         1050.41        9.44e-07            1.56%         2.21e+03     2.61%    2.77e+03          3555   3580     3555   3622
    print(report)
    ind=0
    for i_peak,peak in enumerate(instance_spectrum.list_peaks):
        net_counts_average, uncertainty, total_counts, peak_left_boundary, peak_right_boundary,base_start_left,base_start_right=instance_spectrum.areaTPA(instance_spectrum.ergcal.chn2erg(peak.centroid),width_base=3,base_method=2)
        if net_counts_average>50:
            gamma_match=instance_spectrum.list_gamma_match[i_peak]
            energy=instance_spectrum.ergcal.chn2erg(peak.centroid)
            report_line='{:>2d} \t {:>8d} \t {:>8.2f} \t {:>8s} \t {:>8.2f} \t {:8.2e} \t {:8.2f}% \t {:>10.4e} {:>8.2f}% \t {:>8.2e} \t {:>6d} {:>6d} \t {:>6d} {:>6d}\n'.format(\
                        ind,peak.centroid,energy,gamma_match.nuclide,gamma_match.energy,gamma_match.half_life,gamma_match.intensity*100,net_counts_average, uncertainty*100, total_counts, peak_left_boundary, peak_right_boundary,base_start_left,base_start_right)
            report+=report_line
            print(report_line[:-1])
            ind+=1
    return 0


def export_GV(instance_spectrum,filepath,sample_discription='No sample description was entered.',live_time=10000,total_time=10000):
# def write(list_eff,filepath,increment=0,slope=0.29317,sample_discription='No sample description was entered.',live_time=10000,total_time=10000):
    '''
    Write spectrum files in Spe format. Requires a existing Spe file to dump file

    :param list_eff: --mandatory, counts per channel list of length 8192(for spe format)
    :param filepath: --mandatory, object filepath of spe file
    :param live_time: --optional, live measurement time, 10000 by default
    :param total time: --optional, total measurement time, 10000 by default
    :param sample_discription: --optional, description about sample in string format
    :return: spe formate spectrum in selected filepath
    '''
    spe_open=open(filepath,mode='w',encoding='utf-8')
    spe_open.write('$SPEC_ID:\n')
    spe_open.write(sample_discription+'\n')
    spe_open.write('$SPEC_REM: \n')
    spe_open.write('DET# 1 \n')
    spe_open.write('DETDESC# 1#GB \n')
    spe_open.write('AP# GammaVision Version 6.08 \n')
    spe_open.write('$DATE_MEA: \n')
    time_structure=time.localtime(time.time())
    spe_open.write('%02i/%02i/%04i %02i:%02i:%02i\n'%(time_structure[1],time_structure[2],time_structure[0],time_structure[3],time_structure[4],time_structure[5])) # Needs verification
    spe_open.write('$MEAS_TIM: \n')
    spe_open.write('%5.5s %5.5s\n'%(live_time,total_time))
    spe_open.write('$DATA: \n0 8191\n')    
    i=0
    while True:
        try:
            eff=int(instance_spectrum.list_eff_unsmoothed[i])
            spe_open.write('%8d\n'%(eff))
            i+=1
        except:
            spe_open.write('$ROI:\n')
            spe_open.write('0\n')
            spe_open.write('$PRESETS:\n')
            spe_open.write('None\n0\n0\n')
            spe_open.write('$ENER_FIT:\n')
            spe_open.write('%8.6f %8.6f\n'%(instance_spectrum.ergcal.c0,instance_spectrum.ergcal.c1))
            spe_open.write('$MCA_CAL:\n')
            spe_open.write('3\n')
            spe_open.write('%.6e %.6e %.6e keV\n'%(instance_spectrum.ergcal.c0,instance_spectrum.ergcal.c1,instance_spectrum.ergcal.c2))
            spe_open.write('$SHAPE_CAL:\n')
            spe_open.write('3\n')
            spe_open.write('6.418442E+000 0.000000E+000 0.000000E+000')
            spe_open.close()
            break
    return 0