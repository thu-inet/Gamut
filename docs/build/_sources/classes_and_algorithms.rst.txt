Classes and Algorithms
======================


Following table lists all available classes and their functionalities in GAMUT.

+----------------------------+---------------------------------------------+
|Class                       | Functionality                               |
+============================+=============================================+
|:py:class:`Spectrum`        | store the channel-wise counting data        |
+----------------------------+---------------------------------------------+
|:py:class:`Region`          | represent a ROI on the spectrum which may   |
|                            | contain one or more peaks                   |
+----------------------------+---------------------------------------------+
|:py:class:`Peak`            | represent a peak and store its position,    |
|                            | height, peak width, peak area, etc.         |
+----------------------------+---------------------------------------------+
|:py:class:`Calibration`     | represent the relationship between channel  |
|                            | indexes to variables such as energy         |
+----------------------------+---------------------------------------------+
|:py:class:`Operator`        | the basic class to wrap algorithms          |
+----------------------------+---------------------------------------------+
|:py:class:`Pipe`            | combine multiple operators linearly into one|
|                            | operator                                    |
+----------------------------+---------------------------------------------+
|:py:class:`PipeNet`         | combine multiple operators and pipes        |
|                            | non-linearly into an analysis workflow      |
+----------------------------+---------------------------------------------+
|:py:class:`Node`            | store spectra inside a pipenet              |
+----------------------------+---------------------------------------------+
|:py:class:`Flow`            | wrap operator with its input/output node    |
|                            | information                                 |
+----------------------------+---------------------------------------------+

.. _algotithms:

Following table lists all algorithms implemented in GAMUT.
They are divided into four categories according to their functionalities.
And most algorithms can be replaced by algorithms of the same category easily without any additinal configuration in a defined workflow.


+----------------+---------------------------------------------+
|Functionality   | Algorithms                                  |
+================+=============================================+
|Smoothing       | Savitzky-Golay Fitting Method, Fourier      |
|& Denoising     | Smoothing, Wavelet Smoothing, Whittaker     |
|                | Smoothing Method                            |
+----------------+---------------------------------------------+
|Peak Searching  | Gaussian Factor Method, Numerical           |
|                | Differential Method (Marriscoti Method),    |
|                | Covariance Method, Symmetric Zero Area      |
|                | Method, Double Convolution Method           |
+----------------+---------------------------------------------+
|Baseline        | SNIP method, Self-adaptive SNIP Method,     |
|Stripping       | Iterative Asymmetric Whittaker Stripping    |
+----------------+---------------------------------------------+
|Area Calculation| Total Peak Area Method (Covell,Wasson...),  |
|                | Gaussian Peak Fitting, Deconvolution        |
+----------------+---------------------------------------------+
