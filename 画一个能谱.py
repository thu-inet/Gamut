'''
Author: albertzhang albert.zhangweij@outlook.com
Date: 2024-03-07 00:54:44
Description: 

Copyright (c) 2024 by THU-RSAG, All Rights Reserved. 
'''
import matplotlib.pyplot as plt
import pathlib as pl
import gamspec as gs
import numpy as np

high = gs.Spectrum.from_MCNP(r"E:\Gamut_dev\画一个能谱.c")
# high = gs.Spectrum.from_MCNP(r'E:\Gamut_dev\spectrum_analysis\模拟Eu152能谱\test spectra\Eu152_test_new_10000000.out')
high = gs.Spectrum.from_GammaVision(r'E:\Gamut_dev\spectrum_analysis\标准源谱分析\Eu152_Nie\Eu152 230105.Spe')

high = gs.Spectrum.from_GammaVision(r'D:\User\Spectra\GM20燃耗球\GM20燃耗球\2020-07-21-1800s-023-目标球（反）.Chn')
high = gs.Spectrum.from_GammaVision(r'D:\User\Spectra\GM20燃耗球\GM20燃耗球\2020-07-21-1800s-026-目标球（反）.Chn')

# high = gs.Spectrum.from_GammaVision(r'D:\User\Spectra\GM20燃耗球\GM20燃耗球\2020-07-23-1800s-035-2（反康）.Chn')




plt.figure(figsize=(6, 3), dpi=300)
high.plot(linewidth=0.5, color='gray', chinese_label=True)

# def forward(x):
#     y = np.log10(x)
#     y[x < 1E3] = 3 + (y[x < 1E3] - 3) / 100
#     return y
# def inverse(y):
#     x = 10**y
#     x[y < 3] = 10**((y[y < 3] - 3/2) * 100)
#     return x

# forward = lambda x : x**(1/4)
# inverse = lambda x : x**4

plt.yscale('log')
plt.ylim(1, )
# plt.xlim(0, 1500)
# plt.yscale('function', functions=(forward, inverse))
# plt.ylim(1e3, 1e5)
# plt.legend()
plt.savefig('fig.png')
plt.show()