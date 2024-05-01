![Gamut_logo](https://github.com/Albert-Zhangweijian/Gamut/assets/91001624/eff31c24-9b49-432c-a77a-5d53a189aff3)
# GAMUT
## Introduction

GAMUT (<u>GAM</u>ma <u>UT</u>ils) is an open-source Python framework for gamma spectrum analysis.

Equipped with object-oriented feature, modular design, and pipeline combination capacity, GAMUT is designed to be an accurate and extensible alternative for gamma spectrum analysis, as well as a platform for the implementation of various gamma spectrum analysis algorithms.

GAMUT is now at 0.3.0 version, and more features are under active development. Any trial or secondary development on GAMUT are welcomed! 

## Installation

GAMUT can be installed using Pip with following command:

```shell
pip install gamut-RSAG
```

Or you can clone GAMUT from this repo for easier secondary development:

```shell
git clone https://github.com/thu-inet/Gamut.git
```

## Dependencies

GAMUT relies on following Python packages.

| Package    | Version | Functionality                                  |
| ---------- | ------- | ---------------------------------------------- |
| numpy      |         | Used to realize matrix operation in algorithms |
| pandas     |         | Used to import/export Spectrum objects         |
| matplotlib |         | Used for spectrum visualization                |
| scipy      |         | Used in some algorithms                        |
| pywt       |         | Used in some algorithms                        |

## Quick Start

Here is a MWE (minimum working example) of GAMUT to analysis a gamma spectrum:

```Python
import gamut as gt

# import gamma spectrum
input = gt.Spectrum.from_GammaVision("spectrum.spe")

# define operators(algorithms)
smoother = gt.CentroidSmoother(order=2)

# analysis the spectrum
output = smoother(input)

# plot the spectrum
ax = output.plot()
plt.gcf().savefig("smoothed_spectrum.png")

# export the spectrum
out.export_to_GammaVision("smoothed_spectrum.png")
```

## Classes & Functionalities

Following table lists classes available in GAMUT, 


| Class  |  Functionality|
| -------- | ------------------------------------------------------------ |
| Spectrum    | store the channel-wise counting data                         |
| Region      | represent a ROI on the spectrum which may contain one or more peaks |
| Peak        | represent a peak and store its position, height, peak width, peak area, etc. |
| Calibration | represent the relationship between channel indexes to variables like energy or FWHM |
| Operator | the basic class to wrap algorithms                           |
| Pipe     | combine multiple operators linearly into one operator        |
| PipeNet  | combine multiple operators and pipes non-linearly into an analysis workflow |
| Node     | store spectra inside a Pipenet object                      |
| Flow    | wrap operator with its input/output node information         |

## Reference

[1]  W ZHANG, J LIANG. GAMUT: An Accurate and Extensible Framework for Gamma Spectrum Analysis.

