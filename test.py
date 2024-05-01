'''
Author: albertzhang albert.zhangweij@outlook.com
Date: 2024-01-04 20:15:03
Description: 

Copyright (c) 2024 by THU-RSAG, All Rights Reserved. 
'''
import gamspec
import matplotlib.pyplot as plt

spec = gamspec.Spectrum.from_xml(r'E:\Gamut_dev\模拟谱MixEu000mm_拟合法报告.xml')
spec.plot()
spec.plot_peaks()
plt.show()