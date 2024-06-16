# from .classes.Calibration import *
# from .operator.Operator import *
# from .spectrum.Spectrum import *
# from .classes.Peak import *
# from .classes.Region import *

# from .operator import *

# # branch flow for wavelet smoothing
# wavelet = TranslationInvarianceWaveletSmoother('dmey', 'quadratic-soft', order=3)
# fsmooth = Flow(wavelet, 0, 0, 1, 0)

# # branch flow for peak searching
# savit = SavitzkySmoother(2, 3)
# convol = SecondConvolutionPeakSearcher(10)
# # convol._min_height_ratio = 3
# # convol._min_area_ratio = 3
# # convol._min_height = 20
# # convol._min_area = 100
# searcher = CovarianceSearcher(6, 3, 'uniform')
# searcher._merge_mode = "Old"
# searcher._min_height_ratio = 1
# searcher._min_area_ratio = 1
# searcher._min_height = 0
# searcher._min_area = 0
# fpeak = Flow([savit, convol], 0, 0, 1, 1)

# # merge banch flows
# def combine(spec1, spec2):
#     combined = spec1.copy()
#     combined.regions = deepcopy(spec2.regions)
#     return combined
# comb = Combinor(2, func=combine)
# fcomb = Flow(comb, 1, [0, 1], 2, 0)

# # flow for adaptive SNIP stripper
# strp = OtherOperator.AdaptiveSNIPStripper(baseline=False, high_order=False)
# fstrip = Flow(strp, 2, 0, 3, 0)

# # flow for peak fitting
# fit = RegionPeakFitter(3, equal_width=False, baseline=False)
# ffit = Flow(fit, 3, 0, 4, 0)

# default = PipeNet([fpeak, fsmooth, fcomb, fstrip, ffit])
