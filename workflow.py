# ================================================================
# 模块管理
# ================================================================
import sys
sys.path.append(r"E:\gamspec_dev")
sys.path.append(r"E:\gamspec_dev\gamspec")

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

import gamspec
import gamspec.Spectrum as spe
import gamspec.PeakRegion as pr
import gamspec.operators.Smoother as smr
import gamspec.operators.PeakSearcher as psr
import gamspec.operators.AreaCalculator as acr
import gamspec.operators.OtherOperator as opr
import gamspec.operators.BasicOperator as bop


# ================================================================
# 定义谱分析算子
# ================================================================
sav = smr.SavitzkySmoother(3, 3)
gauss = psr.GaussPeakSearcher(2, 2)
diff = psr.DifferentialPeakSearcher(1, 1, 1, zero_width=4, sigma=0.2)
convol = psr.SecondConvolutionPeakSearcher(sigma=0.1)
strp = opr.SNIPStripper(20)
strp2 = opr.AdaptiveSNIPStripper()
fit = acr.RegionPeakFitter(10, equal_width=False, baseline=False)
minus = bop.Stripper()
calib = opr.FWHMCalibrator()
avg = acr.AverageClassicPeakAreaCalculator(10)
wavelet = smr.TranslationInvarianceWaveletSmoother(mode='soft-quadratic', wavelet='dmey', order=2, step=3)
deconv = acr.Deconvolutioner()
area = acr.AverageClassicPeakAreaCalculator(4)

# ================================================================
# 导入能谱
# ================================================================
# specname = 'MixEu000mm'
# original = gamspec.Spectrum.from_GammaVision(f'D:\\User\\Spectra\\{specname}.Spe')

# specname = 'Co60'
# original = gamspec.Spectrum.from_GammaVision(r'C:\Users\alber\Desktop\CO60 230201 2.5.Spe')

specname = 'Eu152'
original = gamspec.Spectrum.from_GammaVision(r'C:\Users\alber\Desktop\Eu152 230105.Spe')


original.plot()
plt.legend()
plt.savefig(f"{specname}_原始能谱.png")
# plt.show()

# ================================================================
# 定义解谱流程
# ================================================================
smoothed = wavelet(original)  # 平移不变量小波光滑
baseline = strp(smoothed)  # 经典SNIP计算本底谱
striped = minus([smoothed, baseline])  # 剥谱
# searched = diff(striped)  # 数值微分法寻峰
# fitted = fit(searched)  # 拟合计算峰的半高宽

original.plot()
striped.plot()
plt.legend()
plt.savefig(f"{specname}_初次剥谱.png")
# plt.show()

fitted = convol(striped)  # 二次卷积法寻峰
newspec = smoothed.copy()  # 带有本底，但已经光滑
newspec.regions = fitted.regions  # 存取拟合得到的峰数据
baseline2 = strp2(newspec)  # 使用自适应SNIP重新剥谱
striped2 = minus([newspec, baseline2])  # 剥谱


striped2.plot()
striped2.plot_regions()
plt.legend()
plt.savefig(f"{specname}_二次剥谱.png")
plt.show()

# regions_new = []
# witharea = area(striped2)
# for region in witharea.regions:
#     region_new = region.copy()
#     region_new.peaks = []
#     for peak in region.peaks:
#         if peak['area'] > 1000:
#             region.peaks.append(peak)
#     if region_new.peaks:
#         regions_new.append(region_new)
# striped2.regions = regions_new

fitted = fit(striped2)  # 对自适应剥谱的能谱重新拟合
fitted.export_to_xml(f'模拟谱{specname}_拟合法报告.xml')
# classic = avg(striped2)  # 总峰面积法计算峰面积
# calibrated = calib(fitted)  # 根据拟合结果进行半高宽刻度
# deconved = deconv(calibrated)  # 利用刻度结果产生展宽矩阵并反卷积

# ================================================================
# 画图
# ================================================================

# original.plot()
# smoothed.plot()
# baseline.plot()
# striped.plot()
# searched.plot_regions()
# plt.legend()
# plt.savefig(f"{specname}_初步处理.png")
# plt.close()
# # plt.show()

# newspec.plot()
# baseline2.plot()
# striped2.plot()
# striped2.plot_regions()
# plt.legend()
# plt.savefig(f"{specname}_继续处理.png")
# plt.close()
# # plt.show()

# fitted.plot()
# fitted.plot_regions()
# fitted.plot_peaks()
# # deconved.plot(marker='x')
# plt.legend()
# plt.savefig(f"{specname}_拟合结果.png")
# plt.close()

# # deconved.export_to_xml(f'模拟谱{specname}_反卷积法报告.xml')
# fitted.export_to_xml(f'模拟谱{specname}_拟合法报告.xml')
# classic.export_to_xml(f'模拟谱{specname}_总峰面积法报告.xml')

# ================================================================
# 导出分析结果
# ================================================================

df_spectrum = pd.DataFrame()
for spec in [original, smoothed, baseline, striped, fitted, newspec, baseline2, striped2, fitted, ]:
    df_spectrum[spec.label] = spec
df_spectrum.to_excel('模拟谱{specname}_analysis_report.xlsx')