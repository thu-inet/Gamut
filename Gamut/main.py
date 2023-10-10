import numpy as np
import matplotlib.pyplot as plt

import Operator
import Spectrum

import BasicOperator
import OtherOperator

import Smoother
import PeakSearcher
import OtherOperator
import AreaCalculator

# easy test
# cen = Smoother.CentroidSmoother(2)
# strip = OtherOperator.SNIPStripper(4)
# gauss = PeakSearcher.GaussPeakSearcher(2, 0.2)
# wavelet = Smoother.WaveletSmoother('gauss', 'low', 6)
# cal = AreaCalculator.FullPeakAreaCalculator()

# pipe1 = Operator.OperatorPipe(cen, strip, gauss)
# pipe2 = Operator.OperatorPipe(wavelet)
# pipe3 = Operator.OperatorPipe(cal)
# net = Operator.OperatorPipeNet(pipe1, pipe2, pipe3, schemes=[(0, 1), (0, 1), (1, 2)])
# net.plot()

# strip test
simu = Spectrum.SimulatedSpectrum()

sav = Smoother.SavitzkySmoother(3, 4)
strp = OtherOperator.SNIPStripper(5)

lin = BasicOperator.FunctionalOperator(func = lambda x: x)
dualstrp = BasicOperator.Stripper() 

baseline = Operator.Pipe([sav, strp])

schemes = [{'pipe': baseline, 'inp_id': 0, 'inp_order':[0], 'outp_id': 1, 'outp_order': 1},
        {'pipe': lin, 'inp_id': 0, 'inp_order':[0], 'outp_id': 1, 'outp_order': 0},
        {'pipe': Operator.Pipe([dualstrp, sav, sav]), 'inp_id': 1, 'inp_order':[0, 1], 'outp_id': 2, 'outp_order': 0},]
net = Operator.PipeNet(schemes)

net.plot()
net({0:simu})

spec = net.get_node(2)[0]
spec.plot()

simu.plot()
plt.savefig("fig.png")
