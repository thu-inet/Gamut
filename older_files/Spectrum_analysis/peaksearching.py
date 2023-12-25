
class Peak():
    def __init__(self,centroid,left,right,gamma,energy=-1):
        self.centroid=centroid
        self.left=left
        self.right=right
        self.gamma=gamma
        self.energy=energy


def simpleseach(list_eff,prominence,halfwidth):
    list_IF=[1 if list_eff[i]>list_eff[i-1] else 0 for i in range(1,len(list_eff)-1)]
    list_IF.insert(0,0)
    list_IF.append(0)
    i=0
    l=len(list_eff)
    list_peaks=[]
    while 1:
        j=0
        while 1:
            if i+j>=l: 
                return list_peaks
            if list_IF[i+j]==1:
                j+=1
                continue
            else:
                if (j>=halfwidth) & (list_eff[i+j]/max(list_eff[i],1)>=prominence): 
                    peak=Peak(centroid=i+j,left=i,right=min(i+2*j,len(list_eff)),gamma='')
                    list_peaks.append(peak)
                break 
        i+=j+1     


def gauss(list_eff,order=1,factor=1):
    list_PG=[list_eff[i]*list_eff[i+order-2]/list_eff[i-2]/list_eff[i+order]-1 for i in range(2,len(list_eff)-order)]
    list_PG.insert(0,0)
    list_PG.insert(1,0)
    for i in range(order):
        list_PG.append(0)
    list_peaks=[]
    for i in range(len(list_PG)):
        if list_PG[i]>=factor/max(list_eff[i],1)**0.5:
            r=0
            while 1:
                if (list_PG[i+r]+factor/list_eff[i])*(list_PG[i+r+1]+factor/list_eff[i])<=0:
                    break
                else:
                    r+=1
            l=0
            while 1:
                if (list_PG[i-l-1]+factor/list_eff[i])*(list_PG[i-l]+factor/list_eff[i])<=0:
                    break
                else:
                    l-=1
            peak=Peak(centroid=i,left=i-l-1,right=i+r,gamma='')
            list_peaks.append(peak)
    return list_peaks

import matplotlib.pyplot as plt
def drawpeaks(list_erg,list_eff,list_peaks):
    plt.plot(list_erg,list_eff)
    for peak in list_peaks:
        plt.plot(list_erg[peak.left:peak.right],list_eff[peak.left:peak.right],'r-')
    plt.yscale('log')
    plt.show()
    return 0

# if __name__=='__main__':
#     import os
#     import read
#     import smoothing
#     import numpy as np
#     import matplotlib.pyplot as plt

#     filepath=os.path.abspath('.')
#     list_erg,list_eff=read.read(filepath+'\\Spectrum_data\\GM20燃耗球\\2020-07-21-1800s-023-目标球（反）.Chn')
#     list_eff_smthed=smoothing.savgol2(list_eff,5,3)
#     # list_peaks=simpleseach(list_eff_smthed,prominence=1.5,halfwidth=1)
#     list_peaks=gauss(list_eff_smthed,order=2,factor=4)
#     drawpeaks(list_erg,list_eff,list_peaks)



