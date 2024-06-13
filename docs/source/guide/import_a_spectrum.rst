.. _guide_import_a_spectrum:

Import a Spectrum
=================

:class:`Spectrum` is the key class in GAMUT, which is used to store the channel-wise counts of gamma
spectrum. It is essentially a view of :class:`numpy.ndarray`. So it is compatible with all API and
packages based on numpy.

The most authentic way to create a :class:`Spectrum` object is to initiate it with an list-like object.

.. code-block:: python

  # create a sinusoidal spectrum with random noise
  signal = 10 * np.sin(np.linspace(0, 10, 1000))
  signal = signal + np.random.normal(signal.shape)
  spec = gt.Spectrum(signal)

But the more common way is to import it from somewhere else, using corresponding class methods. XML, Pickle, and HDF5 are recommended file formats in GAMUT for that they store all
:class:`Spectrum` object information.

.. code-block:: python

  # import a spectrum form GammVision or other sources
  spec = gt.Spectrum.from_GammaVision("spectrum.spe")
  spec = gt.Spectrum.from_MCNP("mcnp_outp", tally_id=2)
  spec = gt.Spectrum.from_xml("spectrum.xml")
  spec = gt.Spectrum.from_pickle("spectrum.pkl")
  spec = gt.Spectrum.from_hdf5("spectrum.h5")

There are also a class named :class:`SpectrumSimulator`, it is used to generate simulated spectra to
validate algorithm’s performance. Various parameters like length, peaks’ info, and baseline function can be defined. A dictionary named simuspecs has provied some parameter combinations.

.. code-block:: python

  # generate a slightly overlapped doublet
  params = gt.simuspecs['double_slight']
  spec = gt.SimulatedSpectrum(**params)

"""
