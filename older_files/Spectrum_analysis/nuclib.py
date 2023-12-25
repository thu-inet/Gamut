import class_objects
import class_ergcal

def GV2NUIT(nuclide):
    symbol,A=nuclide.split('-')
    if A[-1]=='M':
        A=A[:-1]+'m1'
    if symbol=='J':
        symbol='I'
    return symbol+A

class NYlib():

    def __init__(self,filepath):
        
        with open(filepath) as fileopen:
            filelines=fileopen.readlines()
        ind=0
        while "Nuclide Density" not in filelines[ind]:
            ind+=1
        ind+=8

        self.list=[]
        while "Non-Actinide" not in filelines[ind]:
            line=filelines[ind].split()
            nucid=line[0]
            nuclide='{}'.format(line[1])
            fyield=float(line[-1])
            nuc=class_objects.Nuclide(nuclide,fyield,int(nucid))
            self.list.append(nuc)
            ind+=1

    def match(self,nuclide):
        for nuc in self.list:
            if nuc.nuclide==nuclide:
                return [nuc] 

    def listprint(self,list_nuclide):
        print('Nuclide \t Nucid \t Fyiled\n')
        for nuclide in list_nuclide:
            print('%s    \t %-8.5i \t %-8.5i\n'%(nuclide.nuclide,nuclide.nucid,nuclide.fyield))
        return 0     

class GVlib():

    def __init__(self,filepath):
        self.path=filepath
        with open(filepath,'r',encoding='utf-8') as libopen:
            liblines=libopen.readlines()
        row=6
        self.list=[]
        while True:
            try:
                values=liblines[row].split()
                if values[4]=='Yrs.':
                    half_life=float(values[3])
                elif values[4]=='Days':
                    half_life=float(values[3])/365.25
                elif values[4]=='Hrs.':
                    half_life=float(values[3])/365.25/24
                elif values[4]=='Min.':
                    half_life=float(values[3])/365.25/24/60
                elif values[4]=='Sec.':
                    half_life=float(values[3])/365.25/24/60/60
                else:
                    half_life=0

                gamma=class_objects.Gamma(energy=float(values[0]),
                    intensity=float(values[1][:-1])*0.01,
                    nuclide=values[2],
                    half_life=half_life,
                    NFflags=values[-2],
                    PFflags=values[-1])
                self.list.append(gamma)
                row+=1
            except:
                break

    def listprint(self,list_gamma):
        print('Nuclide \t Energy(keV) \t Intensity(%)\n')
        for gamma in list_gamma:
            print('%s \t %8.5f \t %8.5f\n'%(gamma.nuclide,gamma.energy,gamma.intensity))
        return 0

    def energy_match(self,energy,tolerance=0.005):
        list_diff=sorted(self.list,key=lambda gamma: abs(gamma.energy-energy))
        ind=0
        list_gamma=[]
        while True:
            if abs(energy/list_diff[ind].energy-1)<tolerance:
                list_gamma.append(list_diff[ind])
                ind+=1
            else:
                return list_gamma
    
    def peak_match(self,list_peaks,ergcal,half_life=0,intensity=1E-4,NYlibpath=None,print_FOM=False,tolerance=0.005):
        bias0=0.5
        list_gamma_match=[]
        for peak in list_peaks:
            peak_erg=ergcal.chn2erg(peak.centroid)
            list_candidate=self.energy_match(peak_erg,tolerance=tolerance)
            list_candidate_new=[]
            for gamma in list_candidate:
                if gamma.half_life > half_life and gamma.intensity > intensity and abs(gamma.energy-peak_erg)<5:
                    list_candidate_new.append(gamma)
            list_candidate=list_candidate_new
            if list_candidate==[]:
                list_gamma_match.append(class_objects.Gamma(nuclide='Unknown',energy=0,half_life=0,intensity=0,NFflags='',PFflags=''))
                continue
            
            list_candidate_FOM=[]
            if not NYlibpath:
                for gamma in list_candidate:
                    gamma_FOM=gamma.intensity/gamma.half_life/(abs(gamma.energy-peak_erg)+bias0)
                    list_candidate_FOM.append(gamma_FOM)
            else:
                nylib=NYlib(NYlibpath) # Select by fission yield
                for gamma in list_candidate:
                    gamma_FOM=gamma.intensity/gamma.half_life/(abs(gamma.energy-peak_erg)+bias0)*nylib.match(GV2NUIT(gamma.nuclide))[0].fyield
                    list_candidate_FOM.append(gamma_FOM)
            gamma_match=list_candidate[list_candidate_FOM.index(max(list_candidate_FOM))]
            
            list_gamma_match.append(gamma_match)
            if (print_FOM==True) or (print_FOM==1):
                print('Possible match isotopes for peak at {:8.2f} keV \n'.format(peak_erg))
                print('Nuclide \t Yield \t\t Half life \t Intensity \t Energy bias \t Nuclide FOM \n')
                if NYlibpath:
                    for i_gamma,gamma in enumerate(list_candidate):
                        print('%-7s \t %-.3e \t %-.3e \t %-.3e \t %-.3e \t %-.3e \n'%(gamma.nuclide,nylib.match(GV2NUIT(gamma.nuclide))[0].fyield,gamma.half_life,gamma.intensity,1/(abs(gamma.energy-peak_erg)+bias0),list_candidate_FOM[i_gamma]))
                else:
                    for i_gamma,gamma in enumerate(list_candidate):
                        print('%-7s \t %-.3e \t %-.3e \t %-.3e \t %-.3e \t %-.3e \n'%(gamma.nuclide,0,gamma.half_life,gamma.intensity,1/(abs(gamma.energy-peak_erg)+bias0),list_candidate_FOM[i_gamma]))
        return list_gamma_match

# def reverseGVformat(nuclide):
#     m=1 if nuclide[-1]=='M' else 0
#     symbol,A=nuclide[:len(nuclide)-m].split('-')
#     if symbol=='J':symbol='I'
#     return symbol,int(A),m

# def GVformat(symbol,A,m=0):
#     if m==0:
#         return symbol.capitalize()+'-'+str(A)
#     else:
#         return symbol.capitalize()+'-'+str(A)+'M'

# def NUITformat(symbol,A,m=0):
#     if m==0:
#         return symbol.capitalize()+str(A)
#     else:
#         return symbol.capitalize()+str(A)+'m1'

# def GV2NUIT(nuclide):
#     symbol,A,m=reverseGVformat(nuclide)
#     return NUITformat(symbol,A,m)

    
# if __name__=='__main__':

#     # # Read the nuclib path
#     gvlib=GVlib('D:\\User\\Reports\\Suspect.txt')
#     # print('First 10 gamma lines:\n')
#     # GVlib.listprint(GVlib.list[:10])


#     # Match single gamma line
#     print('Possible match for 661keV:\n')
#     gvlib.listprint(gvlib.energy_match(661))

#     # Read the nuclib path
#     nylib=NYlib('E:\\NUIT_data\\input\\inp_HJ.xml.out')
#     # print('First 10 nuclides:\n')
#     # NYlib.listprint(NYlib.list[:10])

#     # Match nuclides
#     # print('Fyield for Cs-137:\n')
#     # nylib.listprint(nylib.match('Cs137'))

#     # Match goupe gamma line, peaks are example of peaks
#     import peaksearching as ps
#     print('Possible match:\n')
#     list_ergs=[801.58,795.86,604.58,661.37,756.73,724.2,765.8,2185.66,1489.16,696.51,1050.41,873.49,621.93,511.86,427.88,133.51,497.08]
#     list_peaks=[]
#     for erg in list_ergs:
#         peak=ps.Peak(energy=erg,centroid=1,left=1,right=1,gamma=1)
#         list_peaks.append(peak)
#     list_peaks=[ps.Peak(energy=1674.73,centroid=1,left=1,right=1,gamma=1)]
#     # gvlib.listprint(gvlib.peak_match(list_peaks))











