import pywt 
import spectrum as sp
import matplotlib.pyplot as plt
import numpy as np

Cs137=sp.Spectrum()
Cs137.importGV(r'E:\Spectrum_analysis\Spectrum_data\GM20燃耗球\2020-07-21-1800s-023-目标球.Chn')
coeffs=pywt.wavedec(Cs137.list_eff,'sym8','constant',level=4)
print(coeffs)
coeffs_abs=[abs(s) for s in coeffs[-1]]
s=np.median(coeffs_abs)
sigma=s*(2*9)**0.5
for i_cDl in range(1,len(coeffs)):
    for i_cD in range(len(coeffs[i_cDl])):
        # coeffs[i_cDl][i_cD]*=0.5
        detail=coeffs[i_cDl][i_cD]
        if detail > sigma:
            coeffs[i_cDl][i_cD]-=sigma
        elif detail < -sigma:
            coeffs[i_cDl][i_cD]+=sigma
        else:
            pass
print(coeffs)
list_eff_smoothed=pywt.waverec(coeffs,'sym8','constant')
plt.plot(Cs137.list_eff[500:700])
plt.plot(list_eff_smoothed[500:700])
plt.show()