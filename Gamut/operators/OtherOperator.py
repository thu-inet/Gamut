import numpy as np

from ..Spectrum import Spectrum
from ..Operator import Operator
from ..PeakRegion import Calibration


class SNIPStripper(Operator):
    """"
    Classic reverse SNIP baseline stripping algorithm.
    """
    def __init__(self, order: int, label: str = None):
        """
        :param order: time to repeat boardening operation
        :param label: label for the operator
        """
        self._order = order
        if label is None:
            label = "SNIPStripper"
        super().__init__(1, label)

    def __run__(self, spectra: list[Spectrum]) -> Spectrum:
        stripped = spectra[0].copy()
        logloged = np.log(np.log((spectra[0] + 1) ** 0.5 + 1) + 1)
        padded = np.pad(logloged, pad_width=self._order*2, mode='edge')
        for p in range(self._order, 0, -1):
            boardened = np.pad(padded, pad_width=(p, 0), mode='edge')[: -p] + \
                                np.pad(padded, pad_width=(0, p), mode='edge')[p:]
            # padded = np.minimum(padded, boardened / 2)
            padded = np.minimum(padded, boardened / 2)
        sliced = padded[self._order*2: -self._order*2]
        recovered = (np.exp(np.exp(sliced) - 1) - 1) ** 2 - 1
        stripped[:] = recovered
        return stripped


class AdaptiveSNIPStripper(Operator):
    """"
    Classic reverse SNIP baseline stripping algorithm.
    """
    def __init__(self, label: str = None):
        """
        :param order: time to repeat boardening operation
        :param label: label for the operator
        """
        if label is None:
            label = "AdaptiveSNIPStripper"
        super().__init__(1, label)

    def __run__(self, spectra: list[Spectrum]) -> Spectrum:
        stripped = spectra[0].copy()
        logloged = np.log(np.log((spectra[0] + 1) ** 0.5 + 1) + 1)
        for region in stripped.regions:
            order = round(max([peak['stderr']*4 for peak in region.peaks if 'stderr' in peak.keys()]))

            left = min(order, region.left)
            right = min(order, stripped.length-1-region.right)
            val_left = logloged[region.left-left: region.left+1].mean()
            val_right = logloged[region.right: region.right+right+1].mean()
            logloged_region = logloged[region.left-left: region.right+right+1].copy()
            padded_region = np.pad(logloged_region, mode='constant',
                                   constant_values=(val_left, val_right),
                                   pad_width=(order-left, order-right))

            for p in range(order, 0, -1):
                p2 = round(p/2)
                boardened = padded_region[order-p: -order-p] + padded_region[order+p: -order+p+len(padded_region)]
                boardened2 = padded_region[order-p: -order-p] + padded_region[order+p: -order+p+len(padded_region)] \
                    + padded_region[order-p2: -order-p2] + padded_region[order+p2: -order+p2+len(padded_region)]
                padded_region[order: -order] = np.minimum(padded_region[order: -order], np.maximum(boardened/2, boardened2/4))
                # padded_region[order: -order] = boardened2 / 4
            logloged[region.indexes] = padded_region[order: -order].copy()

        recovered = (np.exp(np.exp(logloged) - 1) - 1) ** 2 - 1
        stripped[:] = recovered
        return stripped


class WhittakerSmoother(Operator):
    """
    Whittaker smoother.
    """
    def __init__(self, fit_weight_mode, smooth_weight: float, label: str = None):
        """
        :param order: time to repeat boardening operation
        :param label: label for the operator
        """
        self._fit_weight_mdoe = fit_weight_mode
        self._smooth_weight = smooth_weight
        if label is None:
            label = "WhittakerSmooth"
        super().__init__(1, label)

    def __run__(self, spectra: list[Spectrum], fit_weight: np.ndarray = None, *args, **kargs) -> Spectrum:
        smoothed = spectra[0].copy()
        if self._fit_weight_mdoe == "constant":
            mat_fit_weight = np.diag(np.ones(smoothed.shape))
        elif self._fit_weight_mdoe == "inverse-weighted":
            mat_fit_weight = np.diag(1 / smoothed)
            mat_fit_weight /= mat_fit_weight.mean()
        elif self._fit_weight_mdoe == "self_defined":
            mat_fit_weight = fit_weight
        else:
            raise ValueError("Unknown fit weight mode")
        mat_2nd_diff = np.diag(np.full(smoothed.length, fill_value=2)) + np.diag(np.full(smoothed.length-1, fill_value=-1), k=1) \
            + np.diag(np.full(smoothed.length-1, fill_value=-1), k=-1)
        mat_2nd_diff[0, :3] = [-1, 2, -1]
        mat_2nd_diff[-1, -3:] = [-1, 2, -1]
        smoothed[:] = np.linalg.solve(mat_fit_weight + self._smooth_weight * mat_2nd_diff.T @ mat_2nd_diff, mat_fit_weight @ smoothed)
        return smoothed


class IterativeAsymmetricWhittakerStripper(WhittakerSmoother):

    def __init__(self, smooth_weight: float = 1E3, times: int = 10, label: str = None):
        self._times = times
        if label is None:
            label = "IterativeNonEquilibriumWhittakerStripper"
        super().__init__('self_defined', smooth_weight, label)

    def __run__(self, spectra: list[Spectrum], *args, **kargs) -> Spectrum:
        fit_weight = np.ones(spectra[0].shape[0])
        striped = super().__run__(spectra, np.diag(fit_weight))
        for i in range(self._times):
            fit_weight[striped <= spectra[0]] = 1E-2
            fit_weight[striped > spectra[0]] = 1E-2
            fit_weight[striped > spectra[0]+0.5*spectra[0]**0.5] = 1
            striped = super().__run__(spectra, np.diag(fit_weight))
        return striped


class FWHMCalibrator(Operator):

    def __init__(self, label: str = None):
        if label is None:
            label = "FWHMCalibrator"
        super().__init__(1, label)

    def __run__(self, spectra: list[Spectrum], *args, **kargs) -> Spectrum:
        cal_index, cal_FWHM = [], []
        calibrated = spectra[0].copy()
        for region in calibrated.regions:
            if calibrated.fitness_estimate(region) >= 0.9:
                cal_index.extend([peak['center'] for peak in region.peaks])
                cal_FWHM.extend([peak['stderr']*2.355 for peak in region.peaks])
        calibrated.FWHMcal = Calibration(method='FWHMcal', data=[cal_index, cal_FWHM])
        return calibrated
