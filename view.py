'''
Author: albertzhang albert.zhangweij@outlook.com
Date: 2024-01-09 16:44:28
Description: 

Copyright (c) 2024 by THU-RSAG, All Rights Reserved. 
'''
import gamspec
import gamspec.Spectrum as spe
import matplotlib.pyplot as plt
import pickle

spec = gamspec.Spectrum.from_xml(r'E:\Gamut_dev\模拟谱Eu152_拟合法报告.xml')
origin = gamspec.Spectrum.from_GammaVision(r'C:\Users\alber\Desktop\Eu152 230105.Spe')
origin.regions = spec.regions

with open('test.pkl', 'wb') as f:
    pickle.dump(spec, f)
spec.export_to_excel('test.xlsx')
origin.plot()
origin.plot_peaks()
plt.legend()
plt.show()
