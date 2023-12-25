def write(list_eff,filepath,increment=0,slope=0.29317,sample_discription='No sample description was entered.',live_time=10000,total_time=10000,):
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
            eff=int(list_eff[i])
            spe_open.write('%8d\n'%(eff))
            i+=1
        except:
            spe_open.write('$ROI:\n')
            spe_open.write('0\n')
            spe_open.write('$PRESETS:\n')
            spe_open.write('None\n0\n0\n')
            spe_open.write('$ENER_FIT:\n')
            spe_open.write('%8.6f %8.6f\n'%(increment,slope))
            spe_open.write('$MCA_CAL:\n')
            spe_open.write('3\n')
            spe_open.write('1.605054E+000 2.936760E-001 0.000000E+000 keV\n')
            spe_open.write('$SHAPE_CAL:\n')
            spe_open.write('3\n')
            spe_open.write('6.418442E+000 0.000000E+000 0.000000E+000')
            spe_open.close()
            break
