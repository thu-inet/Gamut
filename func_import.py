import numpy as np
import pathlib
import re

def GammaVision(filepath): 

    if isinstance(filepath, str):
        filepath = pathlib.Path(filepath)
    elif isinstance(filepath, pathlib.Path):
        pass
    else:
        raise TypeError("The filepath should be a string or pathlib.Path() object!")
    
    if filepath.suffix.lower() not in ['.spe', '.chn']:
        raise ImportError("This file cannot be interpreted as Spe format.")
      
    with filepath.open(mode='r', encoding='utf-8') as fileopen:
        filelines = fileopen.readlines()   
        index = 0
    
    # parse spectrum counts
    counts = []
    while not '$DATA' in filelines[index]:
        index += 1
    index += 2
    while '$' not in filelines[index]:
        counts.append( int(filelines[index][:-1]) )
        index += 1

    # parse ROI
    try:
        ROI = []
        while '$ROI' not in filelines[index]:
            index += 1
        index += 2
        while '$' not in filelines[index]:
            m = re.match('([\d]*) ([\d]*)\n', filelines[index])
            ROI.append([int(m.group(1)), int(m.group(2))])
            index += 1
    except:
        ROI = []
    
    # parse energy calibration    
    while '$ENER_FIT' not in filelines[index]:
        index += 1
    index += 1
    m = re.match('([0-9.]*) ([0-9.]*)\n', filelines[index])
    ergcal = {'method': 'linear',
              'c0': float(m.group(2)), 
              'c1': float(m.group(1))}
    
    
    return np.array(counts), ergcal

def MCNP(filepath):

    if isinstance(filepath, str):
        filepath = pathlib.Path(filepath)
    elif isinstance(filepath, pathlib.Path):
        pass
    else:
        raise TypeError("The filepath should be a string or pathlib.Path() object!")
          
    with filepath.open(mode='r', encoding='utf-8') as fileopen:
        filelines = fileopen.readlines()   
        index = 0

    # parse spectrum
    energies, counts = [], []
    while 'tally type 8' not in filelines[index]:
        index += 1
    index += 1
    while 'energy' not in filelines[index]:
        index += 1
    index += 1
    while 'total' not in filelines[index]:
        erg, eff, err = filelines[index].split()
        energies.append( float(erg) )
        counts.append( float(eff) )
        index += 1
    
    # read NPS
    while 'run terminated when' not in filelines[index]:
        index += 1
    m = re.match('\s*run terminated when \s* ([0-9]*) \s* particle histories were done.', filelines[index])
    NPS = int(m.group(1))
    
    return np.array(energies), np.array(counts) * NPS

if __name__ == '__main__':
    
    # eff, cal = GammaVision(filepath=r'E:\Spectrum_analysis\Spectrum_data\GM20燃耗球\2020-07-21-1800s-026-目标球（反）.Spe')
    # erg, eff2 = MCNP(filepath=r'E:\BUMS吸收效应\outputs\outp_inpMIX_100.0000')
    # import matplotlib.pyplot as plt
    # plt.plot(eff)
    # plt.plot(eff2)
    # plt.yscale('log')
    # plt.show()
    