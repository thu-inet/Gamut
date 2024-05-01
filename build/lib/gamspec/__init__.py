'''
Author: albertzhang albert.zhangweij@outlook.com
Date: 2023-10-10 15:57:41
Description: Main entry of Gamut package.

Copyright (c) 2023 by THU-RSAG, All Rights Reserved. 
'''

from .Operator import *
from .Spectrum import *
from .PeakRegion import *

from .operators import BasicOperator as BasicOperator
from .operators import Smoother as Smoother
from .operators import PeakSearcher as PeakSearcher
from .operators import AreaCalculator as AreaCalculator
from .operators import OtherOperator as OtherOperator
# import operators.Smoother as Smoother
# import operators.PeakSearcher as PeakSearcher
# import operators.AreaCalculator as AreaCalculator
# import operators.OtherOperator as OtherOperator
