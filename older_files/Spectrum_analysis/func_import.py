def import_GV(filepath,ergcal): 

    list_erg=[]
    list_eff=[]
    list_eff_unsmoothed=[]

    if not filepath[-3:].lower() in ['spe','chn']:
        raise ValueError('Error-This file cannot be interpreted as Spe or MCNP-ouput format.')
    else:
        pass

    with open(filepath,'r',encoding='utf-8') as fileopen:
        filelines=fileopen.readlines()

    indl=0
    while True:
        line=filelines[indl]
        if '$DATA' in line:
            channels=int(filelines[indl+1].split()[1])
            indl+=2
            break
        else:
            indl+=1

    for i in range(indl,indl+channels):
        line=filelines[i]
        list_eff_unsmoothed.append(int(line[:-1]))
        list_eff.append(int(line[:-1]))
    indl+=channels

    while True:
        line=filelines[indl]
        if '$ENER_FIT' in line:
            indl+=1
            break
        else:
            indl+=1
    
    line=filelines[indl]
    ergcal.c0=float(line.split()[0])
    ergcal.c1=float(line.split()[1])
    list_erg=[ergcal.chn2erg(chn) for chn in range(len(list_eff))]

    return list_erg,list_eff,ergcal

def import_MCNP(filepath):

    list_erg=[]
    list_eff=[]
    list_eff_unsmoothed=[]

    with open(filepath,'r',encoding='utf-8') as fileopen:
        filelines=fileopen.readlines()

    while True:
        line=filelines[indl]
        if 'tally   8' in line:
            indl+=6
            break
        else:
            indl+1

    while True:
        line=filelines[indl]
        try:
            erg,eff,err=line.split()
            list_erg.append(float(erg))
            list_eff_unsmoothed.append(float(eff))
            list_eff.append(float(eff))
            indl+=1
        except:
            break
        
    return list_erg,list_eff