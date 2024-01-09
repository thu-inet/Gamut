'''
Author: albertzhang albert.zhangweij@outlook.com
Date: 2024-01-09 16:44:28
Description: 

Copyright (c) 2024 by THU-RSAG, All Rights Reserved. 
'''
import gamspec
import gamspec.Spectrum as spe
import matplotlib.pyplot as plt

spec = gamspec.Spectrum.from_xml(r'E:\Gamut_dev\模拟谱Eu152_拟合法报告.xml')
origin = gamspec.Spectrum.from_GammaVision(r'C:\Users\alber\Desktop\Eu152 230105.Spe')
origin.regions = spec.regions
origin.plot()
origin.plot_peaks()
plt.legend()
plt.show()
