import numpy as np
import matplotlib.pyplot as plt

def plot_spectrum(spectrum, axe=None, color=None, label=None, linestyle=None):
    
    axe = plt.gca() if axe is None else axe
    color = 'black' if color is None else color
    
    axe.plot(spectrum, linestyle=linestyle, linewidth=1, marker='.', markersize=2, color=color, label=label)  
def plot_peaks(peaks, spectrum, axe=None):
    if axe is None:
        axe = plt.gca()
    for peak in peaks:
        centroid, left, right = peak.centroid, peak.left, peak.right
        axe.fill_between(np.arange(left, right+1), 0, spectrum[left: right+1], color='blue')

def show(ylog=True, axe=None):
    axe = plt.gca() if axe is None else axe
    if ylog == True:
        plt.yscale('log')
    axe.legend()
    plt.show()