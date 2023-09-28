import Operator

import Smoother
import PeakSearcher
import OtherOperator
import AreaCalculator

cen = Smoother.CentroidSmoother(2)
strip = OtherOperator.SNIPstripper(4)
gauss = PeakSearcher.GaussPeakSearcher(2, 0.2)
wavelet = Smoother.WaveletSmoother('gauss', 'low', 6)
cal = AreaCalculator.FullPeakAreaCalculator()

pipe1 = Operator.OperatorPipe(cen, strip, gauss)
pipe2 = Operator.OperatorPipe(wavelet)
pipe3 = Operator.OperatorPipe(cal)
net = Operator.OperatorPipeNet(pipe1, pipe2, pipe3, schemes=[(0, 1), (0, 1), (1, 2)])
net.plot()