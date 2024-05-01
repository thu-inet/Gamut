import numpy as np
from scipy.optimize import fmin, minimize, Bounds
from copy import deepcopy
from time import strftime, ctime
from typing import Literal, Callable

from ..Spectrum import Spectrum
from ..Operator import Operator
from ..PeakRegion import Region

class GenericPeakFitter(Operator):

    def __init__(self, guessfunc: Callable, shapefunc: Callable, renderfunc: Callable,
                 trial: int = 10,
                 label: str | None = None):
        self._trial = trial
        self.guessfunc = guessfunc
        self.shapefunc = shapefunc
        self.renderfunc = renderfunc
        if label is None:
            label = "RegionPeakFitter"
        super().__init__(1, label)

    def __run__(self, spectra: list[Spectrum], *args, **kwargs) -> Spectrum:
        fitted = spectra[0].copy()
        regions = [region for region in fitted.regions if region.npeaks > 0]
        for region in regions:
            print(f"Fitting Region: {region.left}~{region.right}, NPeaks={region.npeaks}, time={strftime(ctime())}")
            best_params, shapes, _ = self._fit(region, fitted, self.guessfunc, self.shapefunc)
            self.renderfunc(region, best_params, shapes)
            print(f"Finish Fitting Region: time={strftime(ctime())}")
        return fitted

    def _fit(self, region: Region, fitted: Spectrum, guessfunc: Callable, shapefunc: Callable):
        best_error = 1E8
        npeaks = region.npeaks
        heights = np.ones(npeaks + 1)  # mutable object to transfer data without return
        fcounts = fitted[region.indexes].copy()
        guesses = guessfunc(region, fitted)
        shapes = shapefunc(region.indexes, guesses)

        def fitfunc(params: np.ndarray, indexes: np.ndarray, counts: np.ndarray, shapefunc: Callable,
                    heights: np.ndarray, shapes: np.ndarray, fcounts: np.ndarray):
            shapes[:] = shapefunc(indexes, params)
            heights[:] = np.abs(np.linalg.lstsq(shapes, counts, rcond=None)[0])
            fcounts[:] = shapes @ heights
            return np.linalg.norm(fcounts - counts, ord=2)
                    
        best_params = guesses.copy()
        for _ in range(self._trial):
            params, error = fmin(fitfunc, x0=guesses,
                                 args=(region.indexes, fitted[region.indexes], shapefunc, heights, shapes, fcounts),
                                 disp=False, full_output=True)[:2]
            if (error <= best_error):
                best_error = error
                guesses = params * np.random.normal(1, 0.01, len(params))
                best_params = params.copy()
            if abs(error/best_error - 1) < 1E-2:
                break
        return best_params, shapes, heights


def tailedgaussian(x: np.ndarray, params: list[float]):
    y = np.exp(-(x-params[0])**2 / (2*params[1]**2))
    y[x < params[0]-params[2]] = np.exp(params[2] * (2*x-2*params[0]+params[2]) / (2*params[1]**2))
    return y
tailedgaussian.nparams = 3

def gaussian(x: np.ndarray, params: list[float]):
    return np.exp(-(x-params[0])**2 / (2*params[1]**2))
gaussian.nparams = 2

def lorentzian(x: np.ndarray, params: list[float]):
    return 1 / ( 1 + (x - params[0]) ** 2 / (params[1] ** 2))
lorentzian.nparams = 2

def linear(x: np.ndarray, params: list[float]):
    return params[0] * x + 1
linear.nparams = 1

def shape_repeater(func: Callable, nshapes: int):
    def repeated(x, params):
        shapes = [func(x, params[i*func.nparams: (i+1)*func.nparams]) for i in range(nshapes)]
        return np.hstack(shapes)
    repeated.nparams = func.nparams * nshapes
    return repeated

def shape_shared_repeater(func: Callable, nshapes: int, shared_param_index: int):
    def shared_repeater(x, params):
        shared_params = params[-1:]
        shapes = []
        for i in range(nshapes):
            i_params = params[i*(func.nparams-1): shared_param_index] + shared_params + params[shared_param_index: (i+1)*(func.nparams-1)]
            shapes.append(func(x, i_params))
        return np.hstack(shapes)
    shared_repeater.nparams = (func.nparams-1) * nshapes + 1
    return shared_repeater

def shape_adder(func: Callable, newfunc):
    def added(x, params):
        return np.hstack((func(x, params[:func.nparams+1]), newfunc(x, params[-newfunc.nparams:])))
    added.nparams = func.nparams + newfunc.nparams
    return added


# class TailedGaussianFitter(GenericPeakFitter):
    
#     def __init__(self, trial: int = 10, label: str | None = None):
#         super().__init__(, self._shape, self._render, trial, label)