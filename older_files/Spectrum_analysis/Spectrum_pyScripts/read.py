import class_ergcal

def read(filepath):
    '''
    Read spectrum files in Spe or MCNP output format.

    :param filepath: --filepath of spe or MCNP-output file
    :return: list_erg and list_eff
    '''
    list_erg=[] 
    list_eff=[]  
    with open(filepath,'r',encoding='utf-8') as fileopen:
        filelines=fileopen.readlines()
    indl=0 # index of file line
    try:
        if filepath[-4:] in ['.spe','.Spe']:
            while True:
                line=filelines[indl]
                if '$ENER_FIT:' in line:
                    slope=float(filelines[indl+1].split()[1])
                    increment=float(filelines[indl+1].split()[0])
                    break
                else:
                    indl+=1
            indl=0
            while True:
                line=filelines[indl]
                indl+=1
                if '0 8191' in line:
                    break
            chn=0 # index of energy channel

            while True:
                try:
                    erg=class_ergcal.energy_calibration(chn,slope=slope,increment=increment)
                    line=filelines[indl]
                    line=line.split()
                    list_eff.append(float(line[0]))
                    list_erg.append(erg)
                    chn+=1
                    indl+=1
                except:
                    del list_erg[-1]
                    del list_eff[-1]
                    break
        else:
            while True:
                line=filelines[indl]
                indl+=1
                if '1tally   8        ' in line:
                    indl+=6
                    break
            while True:
                try:
                    line=filelines[indl]
                    line=line.split()
                    list_erg.append(float(line[0])*1000)
                    list_eff.append(float(line[1]))
                    indl+=1
                except:
                    break
        return list_erg,list_eff
    except:
        print('Error-This file cannot be interpreted as Spe or MCNP-ouput format.')
        return 0



