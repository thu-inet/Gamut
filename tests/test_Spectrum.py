import unittest

from Spectrum import Spectrum
import numpy as np


class TestCaseSepctrum(unittest.TestCase):

    def setUp(self):
        self.spectrum = Spectrum(np.arange(10), 'test')
        self.reference = np.arange(10)

    # 测试取值和切片取值
    def test_getitem(self):
        self.assertEqual(self.spectrum[0], 0)
        self.assertTrue(np.array_equal(self.spectrum[: 5], np.arange(5)))
        self.assertTrue(np.array_equal(self.spectrum[: 6].shape, (6,)))

    # 测试能否像Array一样使用
    def test_attr(self):
        self.assertEqual(self.spectrum.label, 'test')
        self.assertEqual(self.spectrum.mean(), self.reference.mean())
        self.assertEqual(self.spectrum.sum(), self.reference.sum())
        func = lambda x: x ** 2
        self.assertTrue(np.array_equal(func(self.spectrum), func(self.reference)))

    # 测试复制
    def test_copy(self):
        spectrum_copy = self.spectrum.copy()
        self.assertTrue(np.array_equal(self.spectrum, self.reference))
        self.assertEqual(spectrum_copy[0], 0)
        self.assertEqual(spectrum_copy.label, 'test')

if __name__ == '__main__':
    unittest.main()