import matplotlib.pyplot as plt
plt.style.use('bmh')
plt.rcParams['font.family'] = 'Times New Roman'

from matplotlib.pyplot import gca
from numpy import array, arange


class Spectrum:

    def __init__(self, counts:array, label='original', **kargs):
        self.counts = array(counts)
        self.label = label
        self.attrbs = kargs

    def __getattr__(self, attr):
        return getattr(self.counts, attr)

    def __getitem__(self, index):
        return self.counts[index]

    def add(self, spectrum: Spectrum):
        self.counts += spectrum.counts
    
    def sub(self, spectrum):
        self.counts -= spectrum.counts

    def plot(self, *args, ax=None, plot_peaks=False, **kwargs):
        if ax is None:
            ax = gca()
        ax.plot(self.counts, *args, label=self.label, **kwargs)

        if plot_peaks:
            if 'peaks' in self.attrbs.keys():
                for start, end in self.attrbs['peaks']:
                    ax.fill_between(arange(start, end+1), arange(start, end+1)-arange(start, end+1), self.counts[start: end+1])
        if 'ergcal' in self.attrbs.keys():
            ax.set_xlabel('Energy (eV)')
            ax.set_ylabel('Counts')
        else:
            ax.set_xlabel('Channel')
            ax.set_ylabel('Counts')
        ax.set_xlim(0, )
        ax.set_ylim(0, )
        ax.legend()
        return ax

    def copy(self):
        return Spectrum(self.counts, self.label, **self.attrbs)


if __name__ == '__main__':

    import numpy as np
    import unittest

    class TestSpectrum(unittest.TestCase):

        def setUp(self):
            self.spectrum = Spectrum(np.arange(20), 'test', peaks=[(5, 10)])

        # 测试取值和切片取值
        def test_getitem(self):
            self.assertEqual(self.spectrum[0], 0)
            self.assertEqual(self.spectrum[: 5].shape[0], 5)

        # 测试属性和方法
        def test_attr(self):
            self.assertEqual(self.spectrum.label, 'test')
            self.assertEqual(self.spectrum.attrbs['peaks'], [(5, 10)])
            self.assertEqual(self.spectrum.mean(), 9.5)
            self.assertEqual(self.spectrum.sum(), 190)

        # 测试调用画图和其他参数
        def test_call(self):
            ax = self.spectrum.plot()
            self.assertEqual(ax.get_xlabel(), 'Channel')
            self.assertEqual(ax.get_ylabel(), 'Counts')

        # 测试复制
        def test_copy(self):
            spectrum_copy = self.spectrum.copy()
            self.assertEqual(spectrum_copy[0], 0)
            self.assertEqual(spectrum_copy.label, 'test')
            self.assertEqual(spectrum_copy.attrbs['peaks'], [(5, 10)])

        # 测试深拷贝
        def test_copy_change(self):
            spectrum_copy = self.spectrum.copy()
            spectrum_copy.counts[0] = 1
            self.assertEqual(self.spectrum[0], 0)
            spectrum_copy.label = 'test2'
            self.assertEqual(self.spectrum.label, 'test')
        
        # 其他测试
        def test_others(self):
            spectrum2 = Spectrum(np.ones(20), 'test2')
            # self.assertEqual(self.spectrum - spectrum2, np.arange(20) - np.ones(20))

    unittest.main()
