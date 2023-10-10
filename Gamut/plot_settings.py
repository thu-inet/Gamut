import matplotlib as mpl
import matplotlib.pyplot as plt


_config = {'figure.dpi': 200,
           'figure.figsize': (6, 4),
           'font.size': 8,
           'lines.linewidth': 1,
           'lines.markersize': 4,
           'lines.markeredgewidth': 0.5,
           'font.family': 'Times New Roman'}
_colors = plt.get_cmap('tab20').colors
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
