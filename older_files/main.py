
import matplotlib.pyplot as plt

spectrum = Spectrum(simu_spectrum.simu_spectrum(length=1000))

centroid = CentroidSmoother(10)
savitzky = SavitzkySmoother(3, 5)
fourier = FourierSmoother('low', 0.1)
pipe_censav = OperatorPipe(centroid, savitzky)
pipe_foufou = OperatorPipe(fourier, fourier)
pipe_censavfou = OperatorPipe(centroid, savitzky, fourier)
voidpipe = OperatorPipe()  # test void pipe

spe_cen = centroid(spectrum)
spe_sav = savitzky(spectrum)
spe_fou = fourier(spectrum)
spe_pipe = pipe_censav(spectrum)

# ax = spectrum.plot()
# ax_cen = spe_cen.plot()
# ax_sav = spe_sav.plot()
# ax_fou = spe_fou.plot()
# ax_pipe = spe_pipe.plot()
# plt.legend(fontsize=6)
# plt.show()

# test pipenet
# schemes = [(0, 1), (0, 2), (2, 1),
#            (1, 8), (1, 7), (2, 5),
#            (5, 7), (5, 8), (7, 8)]
# pipenet = OperatorPipeNet(pipe_censavfou, pipe_foufou, pipe_censav,
#                           pipe_censav, pipe_foufou, pipe_censav,
#                           pipe_censav, pipe_foufou, pipe_censav, schemes=schemes)
# pipenet(spectrum)
# pipenet.plot()
# for spectrums in pipenet._pipedspectrums:
#     for spectrum in spectrums:
#         ax = spectrum.plot()
# plt.legend()
# plt.show()

# test MBC
MBCcentroid = MBCTransformer(CentroidSmoother)(10)
spe_mbccen = MBCcentroid(spectrum)
ax_cen = spe_cen.plot()
ax_mbccen = spe_mbccen.plot()
ax_mbccen.set_yscale('linear')
plt.legend()
plt.show()
