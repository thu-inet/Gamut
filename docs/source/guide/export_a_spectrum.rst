.. _guide_export_a_spectrum:

Export a Spectrum
=================

A Spectrum object can be exported as several file formats, but some information may need to be
added manually, such as the measurement time when exported into a GammaVision output file.
Yet XML, Pickle and HDF5 store exactly all data of a Spectrum object and are recommended.

.. code-block:: python

  spec.export_to_GammaVision('analyzed.spe')
  spec.export_to_Pickle('analyzed.pkl')

Furthermore, Spectrum has a method called export_to_pandas(), which only exports the peaks’
data into a pandas.DataFrame. It is convinient to generate the final analysis results.

.. code-block:: python

  # save all peaks' data into an Excel file
  spec.export_to_pandas().to_excel('results.xlsx')
