import numpy as np

from copy import deepcopy
from typing import Optional

from .Peak import Peak
from ..utils import gaussian

class Region:
    """
    Class representing a region of interest in a spectrum.
    It requires the left and right indexes of the region, and can contain multiple peaks.
    Parameters like slope or offset are optional, depending if the baseline has been fitted.
    Other trivial parameters are optional, and are stored in the attrs dictionary.
    """
    def __init__(self, left: int,right: int,
                 peaks: Optional[list[Peak]] = None,
                 slope: Optional[float] = None,
                 offset: Optional[float] = None, **kwargs) -> None:
        """
        Initialize the region object with the given left and right indexes and optional parameters.
        
        :param left: The left index of the region
        :param right: The right index of the region
        :param peaks: The peaks in the region
        :param slope: The slope of the baseline
        :param offset: The offset of the baseline
        """
        self._left = left
        self._right = right
        self._peaks = peaks if peaks is not None else []
        self._slope = slope
        self._offset = offset
        self.attrs = kwargs

    @property
    def left(self) -> int:
        """
        The left index of the region.
        """
        return self._left
    
    @left.setter
    def left(self, left: int) -> None:
        self._left = left
        
    @property
    def right(self) -> int:
        """
        The right index of the region.
        """
        return self._right
    
    @right.setter
    def right(self, right: int) -> None:
        self._right = right
    
    @property
    def peaks(self) -> list[Peak]:
        """
        The list of peaks in the region.
        """
        return self._peaks[:]

    @peaks.setter
    def peaks(self, peaks: list[Peak]) -> None:
        self._peaks = peaks[:]    

    @property
    def npeaks(self) -> int:
        """
        Number of peaks in the region.
        """
        return len(self.peaks)

    @property
    def length(self) -> int:
        """
        Length of the region.
        """
        return self.right - self.left + 1

    @property
    def indexes(self) -> np.ndarray:
        """
        Indexes of the region.
        """
        return np.arange(self.left, self.right+1)
    
    @property
    def slope(self) -> Optional[float]:
        """
        The slope of the baseline. Available only if the baseline has been fitted.
        """
        return self._slope
    
    @slope.setter
    def slope(self, slope: float) -> None:
        self._slope = slope
        
    @property
    def offset(self) -> Optional[float]:
        """
        The offset of the baseline. Available only if the baseline has been fitted.
        """
        return self._offset
    
    @offset.setter
    def offset(self, offset: float) -> None:
        self._offset = offset

    def __repr__(self) -> str:
        return f'Region({self.left}~{self.right})[{self.npeaks}]'    

    def copy(self) -> 'Region':
        """
        Return a deep copy of the region
        """
        return deepcopy(self)
    
    def shift_by_index(self, shift: int) -> None:
        """
        Shift the region position by a given index, this will also shift the peaks.
        
        :param shift: int, the shift value
        """
        self._left += shift
        self._right += shift
        for peak in self.peaks:
            peak.shift_by_index(shift)

    def fit_peak(self, peak: Peak, indexes: Optional[np.ndarray] = None) -> np.ndarray:
        """
        Return the fitted spectrum segment using the peak parameters.
        
        :param peak: The peak to fit
        :param indexes: The indexes of the region

        :return: The fitted spectrum segment
        """
        indexes = indexes if indexes is not None else self.indexes
        if (peak.center is not None) and (peak.height is not None) and (peak.stderr is not None):
            fitted_peak = gaussian(indexes, peak.center, peak.stderr) * peak.height
        else:
            raise ValueError('Peak parameters are not complete.')
        return fitted_peak

    def fit_baseline(self, indexes: Optional[np.ndarray] = None) -> np.ndarray:
        """
        Return the fitted baseline using the baseline parameters.
        
        :param indexes: The indexes of the region
        
        :return: The fitted baseline
        """
        indexes = indexes if indexes is not None else self.indexes
        if self.offset is not None and self.slope is not None:
            fitted_baseline = indexes * self.slope + self.offset
        else:
            raise ValueError('Baseline parameters are not complete.')    
        return fitted_baseline

