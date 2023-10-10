import unittest

from Operator import *
from Spectrum import *
import numpy as np


class TestCaseOperator(unittest.TestCase):

    def setUp(self) -> None:
        # Single Input Operator Initialization
        self.linear = FunctionalOperator(func=lambda x: x, label='linear')
        self.modulo = FunctionalOperator(func=lambda x: x % 4, label='modulo')
        self.divide = FunctionalOperator(func=lambda x: x / 4, label='divide')

        # Multiple Input Operator Initialization
        self.multiply = Multiplier()

        # Test Spectrum Initialization
        self.spectrum = [Spectrum(counts=np.arange(1, 5), label='spe1')]
        self.spectra = [Spectrum(counts=np.arange(1, 5), label='spe1'), Spectrum(counts=np.arange(2, 6), label='spe2')]

    def test_linear_output(self):
        self.assertTrue(np.array_equal(self.linear(self.spectrum), np.arange(1, 5)))

    def test_modulo_output(self):
        self.assertTrue(np.array_equal(self.modulo(self.spectrum), np.arange(1, 5) % 4))

    def test_divide_output(self):
        self.assertTrue(np.array_equal(self.divide(self.spectrum), (np.arange(1, 5) / 4)))

    def test_multiply_output(self):
        self.assertTrue(np.array_equal(self.multiply(self.spectra), np.arange(1, 5)*np.arange(2, 6)))

    def test_form_pipe(self):
        Pipe([self.linear, self.modulo, self.divide])
        self.assertRaises(AttributeError, Pipe, [self.linear, self.modulo, self.multiply])
        Pipe([self.multiply, self.linear, self.modulo])

if __name__ == '__main__':
    unittest.main()
