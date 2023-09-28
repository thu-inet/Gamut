'''
 # @ Author: Weijian ZHANG
 # @ Create Time: 2023-07-19 16:23:15
 # @ Modified by: Weijian ZHANG
 # @ Modified time: 2023-07-19 16:23:29
 # @ Description:
 '''

from Spectrum import Spectrum
from Operator import Operator
from abc import ABC, abstractmethod


# class CountingAreaCalculator(Operator):
    
#     def __init__(self, label='area'):
#         super().__init__(label)
#         pass
    
#     @ abstractmethod
#     def __run__(self, spectrum):
#         pass

class FullPeakAreaCalculator(Operator):
    
    def __init__(self, label='-FPAreaCalculator[Convell]'):
        super().__init__(label)
        pass

    def __run__(self, spectrum):
        pass