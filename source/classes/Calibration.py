import numpy as np

from typing import Tuple
from scipy.optimize import curve_fit

from ..globals import CalibrationMethod, StructuredFunction, calibration_functions


class Calibration:
    """
    Relationship between spectrum channels and other quantities.
    It uses a funtion and a set of parameters to describe the relationship.
    When the parameters are unknown, they can be fitted using a set of points.
    """

    def __init__(self, method: CalibrationMethod,
                 params: list[float]):
        """
        Initialize the calibration object with the given method and parameters.
        
        :param method: The calibration method to use
        :param params: The parameters of the calibration function
        """
        self._method = method
        self._params = params
        self._points = [] # if directly defined, then no data is needed

    @property
    def method(self) -> CalibrationMethod:
        """
        Alias of the curve fitting function.
        """
        return self._method

    @method.setter
    def method(self, method: CalibrationMethod) -> None:
        self._method = method

    @property
    def params(self) -> list[float]:
        """
        Parameters used in the curve fitting function.
        """
        return self._params

    @params.setter
    def params(self, params: list[float]) -> None:
        self._params = params

    @property
    def function(self) -> StructuredFunction:
        """
        Structured function used in the curve fitting.
        """
        return calibration_functions[self.method]

    @property
    def points(self) -> list[Tuple[float, float]]:
        """
        List of points(channel, quantity) for curve fitting.
        """
        return self._points.copy()

    @points.setter
    def points(self, points: list[Tuple[float, float]]) -> None:
        self._points = points.copy()
        self._params = self._fitparams(self.method, self.points)

    @property
    def points_indexes(self) -> list[float]:
        """
        Indexes of the points used for calibration.
        """
        return [item[0] for item in self.data]

    @property
    def points_values(self) -> list[float]:
        """
        Values of the points used for calibration.
        """
        return [item[0] for item in self.data]
    
    def add_point(self, point: Tuple[float, float]) -> None:
        """
        Add a point to the points list.
        """
        self._points.append(point)
        self._params = self._fitparams(self.method, self.points)

    def __call__(self, inp: float | np.ndarray) -> float | np.ndarray:
        return self.function.function(inp, self.params)

    @classmethod
    def from_points(cls, method: CalibrationMethod, points: list[Tuple[float, float]]) -> 'Calibration':
        """
        Initialize the calibration object with the given method and points.
        The function parameters are fitted using these points.
        
        :param method: The calibration method to use
        :param points: The points used for curve fitting
        """
        params = cls._fitparams(method, points)
        return cls(method, params)

    @staticmethod
    def _fitparams(method: CalibrationMethod, data: list[Tuple[float, float]]) -> Tuple[list[float], float]:
        """
        Fit the function parameters using the given points.

        :param method: The calibration method to use
        :param data: The points used for curve fitting
        """
        f = calibration_functions[method]
        indexes = np.array([item[0] for item in data])
        values = np.array([item[1] for item in data])
        params, _ = curve_fit(f, xdata=indexes, ydata=values)
        return params
