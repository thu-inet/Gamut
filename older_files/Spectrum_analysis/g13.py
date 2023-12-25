from spectrum import *

def readg13(filepath=r'C:\users\Alber\Desktop\g13'):
    with open(filepath,'r') as fileopen:
        filelines=fileopen.readlines()

    spectrum=[[],[],[]]
    i=2
    while True:
        try:
            info=filelines[i].split()
            spectrum[0].append(float(info[0])*1000)
            spectrum[1].append(float(info[1]))
            spectrum[2].append(float(info[2]))
            i+=1
        except:
            break
    
    return spectrum

if __name__ == '__main__':

    filepath=r'E:\WeChat Files\wxid_vebikycea99422\FileStorage\File\2022-11\g7'
    list_erg,list_eff,list_eff_simu=readg13(filepath=filepath)
    g13=Spectrum(list_erg=list_erg,list_eff=list_eff,list_eff_unsmoothed=list_eff)
    g13.ergcal.recalibrate(list_chn=[3657,15997],list_erg=[list_erg[3657],list_erg[15997]],L=15998)
    
    g13.smooth_fourier(factor=0.3)
    g13.search_simple()
    # g13.add_peak(peak_left_boundary= 3656,peak_right_boundary=3662)
    # Cs137.singlematch(1384,NYlibpath='E:\\NUIT_data\\input\\inp_test.xml.out')
    # g13.match(NYlibpath='E:\\NUIT_data\\input\\inp_HJ.xml.out')
    g13.draw_spectrum()
    # g13.report()


    

