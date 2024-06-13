from typing import Optional
from copy import deepcopy

class Peak:
    """
    Class representing a single peak in a spectrum.
    The only required parameter is the location/index of the peak.
    Parameters like center or height are optional, depending if the peak has been fitted. 
    Other trivial parameters are stored in the attrs dictionary.
    """
    def __init__(self, location: int,
                 center: Optional[float] = None,
                 height: Optional[float] = None,
                 stderr: Optional[float] = None, **kwargs) -> None:
        """
        Initialize the peak object with the given location and optional parameters.
        
        
        :param location: The index of the peak in the spectrum
        :param center: The center of the peak
        :param height: The height of the peak
        """
        
        self._location = location
        self._center = center
        self._height = height
        self._stderr = stderr
        self.attrs = kwargs

    @property
    def location(self) -> int:
        """
        The center index of the peak.
        """
        return self._location

    @location.setter
    def location(self, location: int) -> None:
        self._location = location

    @property
    def center(self) -> Optional[float]:
        """
        The precise center of the peak.
        """
        return self._center
    
    @center.setter
    def center(self, center: float) -> None:
        self._center = center
    
    @property
    def height(self) -> Optional[float]:
        """
        The height of the peak.
        """
        return self._height
    
    @height.setter
    def height(self, height: float) -> None:
        self._height = height
    
    @property
    def stderr(self) -> Optional[float]:
        """
        The standard error of the gaussian shape of the peak.
        """
        return self._stderr
    
    @stderr.setter
    def stderr(self, stderr: float) -> None:
        self._stderr = stderr
        
    def copy(self) -> 'Peak':
        """
        Return a deep copy of the peak.
        """
        return deepcopy(self)

    def shift_by_index(self, shift: int) -> None:
        """
        Shift the peak position by a given index
        
        :param shift: Shift value on the index, positive indicates the right direction.
        """
        self._location += shift
        self._center += shift if self._center is not None else None