import numpy as np
from copy import deepcopy
from numpy.typing import NDArray
from typing import Callable

# from ..operator import Operator
from ..classes import Region, Calibration
from . import *

class Spectrum(np.ndarray):
    """
    Spetrum is the most basic object of GAMUT. It is used to store data of a gamma spectrum.
    It is a view of 1-d numpy.ndarray.
    """
    def __init__(self,
                 counts: list[float] | NDArray[np.float64],
                 label: str | None = None,
                 regions: list[Region] | None = None,
                 attrs: dict | None = None,
                 ergcal: Calibration | None = None,
                 FWHMcal: Calibration | None = None):
        """
        Initialize the spectrum object with the given counts and optional parameters.
        
        :param counts: The channel-wise counts of the spectrum
        :param label: Optional, the name or label of the spectrum
        :param regions: Optional, the list of ROIs of the spectrum
        :param attrs: Optional, the dictionary of trivial attributes of the spectrum
        :param ergcal: Optional, the energy calibration of the spectrum
        :param FWHMcal: Optional, the FWHM calibration of the spectrum
        """
        self._label = label if label is not None else 'spectrum'
        self._attrs = attrs if attrs is not None else {}
        self._regions = regions if regions is not None else []
        self._ergcal = ergcal  # let the calibration be None if the spectrum is not calibrated
        self._FWHMcal = FWHMcal

    @property
    def label(self) -> str:
        """
        The name or label of the spectrum.
        """
        return self._label
    
    @label.setter
    def label(self, label: str) -> None:
        self._label = label
        
    @property
    def regions(self) -> list[Region]:
        """
        The list of ROIs in the spectrum.
        """
        return self._regions
    
    @regions.setter
    def regions(self, regions: list[Region]) -> None:
        self._regions = regions
        
    @property
    def attrs(self) -> dict:
        """
        The dictionary of trivial attributes of the spectrum.
        """
        return self._attrs
    
    @attrs.setter
    def attrs(self, attrs: dict) -> None:
        self._attrs = attrs
    
    @property
    def ergcal(self) -> Calibration:
        """
        The energy calibration of the spectrum.
        """
        return self._ergcal
    
    @ergcal.setter
    def ergcal(self, ergcal: Calibration) -> None:
        self._ergcal = ergcal
    
    @property
    def FWHMcal(self) -> Calibration:
        """
        The FWHM calibration of the spectrum.
        """
        return self._FWHMcal

    @FWHMcal.setter
    def FWHMcal(self, FWHMcal: Calibration) -> None:
        self._FWHMcal = FWHMcal

    @property
    def length(self) -> int:
        """
        The length of the spectrum.
        """
        return self.shape[0]

    @property
    def indexes(self) -> np.ndarray:
        """
        The indexes of the spectrum.
        """
        return np.arange(self.length)

    @property
    def energies(self) -> np.ndarray:
        """
        The energies of the spectrum. Available only when the spectrum is calibrated.
        """
        if self.ergcal is None:
            raise ValueError("The spectrum is not calibrated.")
        return self.ergcal(self.indexes)

    @property
    def counts(self) -> np.ndarray:
        """
        The channel-wise counts of the spectrum.
        """
        return np.asarray(self)

    @property
    def nregions(self) -> int:
        """
        The number of regions in the spectrum.
        """
        return len(self.regions)

    @property
    def npeaks(self) -> int:
        """
        The number of peaks in the spectrum.
        """
        return sum([len(region.peaks) for region in self.regions])

    @property
    def peaks(self):
        """
        The list of peaks in the spectrum. It retuns a copy of the peaks.
        """
        return [peak for region in self.regions for peak in region.peaks].copy()

    def __rshift__(self, operator: Callable):
        """
        Overload the right shift operator as the pipeline operation.
        
        :param operator: The operator to be applied to the spectrum.
        """
        return operator(self) # overload the right shift operator as the pipeline operation

    def __str__(self):
        return f"Spectrum[{self.label}]"
    
    def from_GammaVision(self, filename: str):
        """
        Import the spectrum from a GammaVision file.
        
        :param filename: The path to the GammaVision file.
        """
        return SpectrumImporter.from_GammaVision(filename)
    
    def from_excel(self, filename: str, counts_column: str, energy_column: str | None = None):
        """
        Import the spectrum from columns of excel files.
        
        :param filename: The path to the excel file.
        :param counts_column: The name of the column containing the counts.
        :param energy_column: The name of the column containing the energies.
        """
        return SpectrumImporter.from_excel(filename, counts_column, energy_column)

    # ==========================================================================
    # special methods to make the subclassing of numpy.ndarray work properly
    # ==========================================================================
    def __new__(cls, counts, *args, **kargs):
        self = np.asarray(counts, dtype=np.float64).view(cls)
        return self

    def __array__wrap__(self, array, context=None):
        if array.ndim == 0:
            return float(array)
        else:
             return np.ndarray.__array_wrap__(self, array, context)

    def __array_finalize__(self, obj):
        self.label = deepcopy(getattr(obj, 'label', 'spectrum'))
        self.ergcal = deepcopy(getattr(obj, 'ergcal', None))
        self.FWHMcal = deepcopy(getattr(obj, 'FWHMcal', None))
        self.attrs = deepcopy(getattr(obj, 'attrs', {}))
        self.regions = deepcopy(getattr(obj, 'regions', []))

    def __reduce__(self):
        pickled_state = super().__reduce__()
        if isinstance(pickled_state[2], tuple):
            new_state = pickled_state[2] + (self.label, self.regions, self.ergcal, self.FWHMcal, self.attrs)
            return (pickled_state[0], pickled_state[1], new_state)
        else:
            raise ValueError("pickled_state[2] is not a tuple(Unsolved problem)")

    def __setstate__(self, state):
        self.label = state[-5]
        self.regions = state[-4]
        self.ergcal = state[-3]
        self.FWHMcal = state[-2]
        self.attrs = state[-1]
        super().__setstate__(state[0:-5])

    def __getitem__(self, key):

        sliced = super().__getitem__(key)
        if not isinstance(key, slice):
            return sliced

        start = key.start if key.start is not None else 0
        stop = key.stop if key.stop is not None else self.length

        if self.ergcal.by == 'data':
            data = self.ergcal.data.copy()
            data[:, 0] -= start
            ergcal = Calibration(method=self.ergcal.method, data=data)
        else:
            params = self.ergcal.params.copy()
            params[0] += self.ergcal(start)
            ergcal = Calibration(method=self.ergcal.method, params=params)

        regions = []
        for region in self.regions:
            if (region.right <= start) or (region.left >= stop):
                continue
            else:
                region_copy = region.copy()
                region_copy.left = max(region_copy.left-start, 0)
                region_copy.right = min(region_copy.right-start, stop-start)
                for peak in region_copy.peaks:
                    peak['location'] = min(max(peak['location']-start, 0), stop-start)
                    if 'center' in peak.keys():
                        peak['center'] = min(max(peak['center']-start, 0), stop-start)
                regions.append(region_copy)

        return Spectrum(sliced, label=self.label, regions=regions, ergcal=ergcal, FWHMcal=self.FWHMcal, attrs=self.attrs)
    

