import old.read as read
import areaTPA

filepath=r'E:\\Spectrum_Analysis\\Data\\'
list_erg,list_eff=read.read(filepath+'\\Cs-137.spe')
list_ergb,list_effb=read.read(filepath+'\\background.spe')

area,error,total,left,right,basel,baser=areaTPA.areaTPA(list_erg,list_eff,peak=2255,width_base=20,avg_number=10)
print('TPA Peakarea：',area)
print('TPA Uncertainty：',error)
print('TPA Total counts：',total)
print('TPA Peak_boundary_left：',left)
print('TPA Peak_boundary_right：',right)
# areab,errorb,totalb,leftb,rightb,d,d2=areaTPA.areaTPA(list_erg,list_effb,peak=2657)
# print('TPA Peakarea：',areab)
# print('TPA Uncertainty：',errorb)
# print('TPA Total counts：',totalb)
# print('TPA Peak_boundary_left：',leftb)
# print('TPA Peak_boundary_right：',rightb)

l=left
r=right
def deduc(list_eff,list_effb,l,r):
    l2=len(list_eff)
    list_eff_smoothed=[]
    list_effb_smoothed=[]
    for i in range(l2):
        # list_eff_smoothed.append((2*list_eff[i]+list_eff[max(i-1,0)]+list_eff[min(i+1,l-1)])/4) # 3-points smoothing method  
        # list_eff_smoothed.append((4*list_eff[i]+list_eff[max(i-1,0)]+list_eff[min(i+1,l-1)])/6) # 3-points Simpson smoothing method 
        list_eff_smoothed.append((-3*list_eff[max(i-2,0)]+12*list_eff[max(i-1,0)]+17*list_eff[i]+12*list_eff[min(i+1,l2-1)]-3*list_eff[min(i+2,l2-1)])/35) # # 5-points smoothing method 
        list_effb_smoothed.append((-3*list_effb[max(i-2,0)]+12*list_effb[max(i-1,0)]+17*list_effb[i]+12*list_effb[min(i+1,l2-1)]-3*list_effb[min(i+2,l2-1)])/35) # # 5-points smoothing method 
    list_eff=list_eff_smoothed
    list_effb=list_effb_smoothed
    total=sum(list_eff[l:r+1])
    back=sum(list_effb[l:r+1])*88562/62139
    error=(total+back*88562/62139)**0.5
    return total,back,total-back,error
print(deduc(list_eff,list_effb,l,r))