
import matplotlib as mpl
import matplotlib.pyplot as plt

from typing import Literal, TypeAlias
from collections import namedtuple

# function with fixed number of parameters
# convention: the first parameter is the offset
StructuredFunction = namedtuple('CalibrationFunction', ['name', 'nparams', 'function'])
CalibrationMethod: TypeAlias = Literal['linear', 'quadratic', 'FWHM']
calibration_functions = {
    'linear': StructuredFunction('linear', 2, lambda x, params: x * params[1] + params[0]),
    'quadratic': StructuredFunction('quadratic', 3, lambda x, params: x * params[2] ** 2 + params[1] * x + params[0]),
    'FWHM': StructuredFunction('FWHM', 3, lambda x, params: params[0] + params[1] * (x + params[2] * x ** 2) ** 0.5)}



# configure the matplotlib
_config = {'figure.dpi': 200,
           'figure.figsize': (6, 4),
           'font.size': 10,
           'lines.linewidth': 1,
           'lines.markersize': 4,
           'lines.markeredgewidth': 0.5,
           'font.family': ['Times New Roman', 'FangSong']}
_colors = plt.get_cmap('Set3').colors
_markers = ['o', 'v', '^', 's', 'p', '*', 'h', 'H', 'D', 'd']

mpl.rcdefaults()
plt.style.use('bmh')
mpl.rcParams.update(_config)

plt.rcParams['axes.labelsize'] = plt.rcParams['font.size']
plt.rcParams['axes.titlesize'] = 1.5*plt.rcParams['font.size']
plt.rcParams['legend.fontsize'] = plt.rcParams['font.size']
plt.rcParams['xtick.labelsize'] = plt.rcParams['font.size']
plt.rcParams['ytick.labelsize'] = plt.rcParams['font.size']


def markers(i: int):
    return _markers[i % len(_markers)]


def colors(i: int):
    return _colors[i % len(_colors)]




