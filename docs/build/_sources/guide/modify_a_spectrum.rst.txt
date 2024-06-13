Modify a Spectrum
==================

Aside from the above-mentioned usage, sometimes users may want to modify the spectra or other
objects manually. It can be done easily as most classes are inherited from common classes in
Python. But it should be careful in such modification.

.. code-block:: python

  # add a ROI
  region = gt.Region(left=40, right=60)
  spec.regions.append(region)

  # remove the third peak in a ROI
  region.remove(2)# slice the spectrum

  # the index of ROIs and peaks will be changed automatically
  sliced = spec[100: 200]
