import unittest
import numpy as np
from Gamut.Spectrum import *


class TestSpectrum(unittest.TestCase):

    def setUp(self):
        self.spectrum = Spectrum(np.array([1, 2, 3]))

    def test_init(self):
        self.assertEqual(Spectrum(np.array([1, 2, 3])).counts, np.array([1, 2, 3]))

    def tearDown(self):
        return super().tearDown()

if __name__ == '__main__':
    # unittest.main()
    suite = unittest.TestLoader().loadTestsFromTestCase(TestSpectrum)
    unittest.TextTestRunner(verbosity=2).run(suite)
