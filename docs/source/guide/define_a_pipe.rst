
.. currentmodule:: gamut.operator
.. currentmodule:: gamut.spectrum

Define a Pipe
=============


When multiple algorithms should be executed in sequence, they can be concatenated into a :py:class:`Pipe`
object. A Pipe object is just like an :py:class:`Operator` object. When it is called, it executes its first
internal :py:class:`Operator` object, and then passes its output spectrum to the next, and so on,
until all :py:class:`Operator` objects are called once.

The :py:attr:`inp_num` of all subsequent operators must equals to 1, and the inp_num of the :py:class:`Pipe` equals
to that of its first :py:class:`Operator` object.


.. code-block:: python

  # define two peak searching operators
  covol = gt.DoubelConvolutionPeakSearcher(hwidth=10)
  covar = gt.CovariancePeakSearcher(hwidth=10, FWHM=4, weight_mode='inverse')
  covar.merge_mode = 'Old'
  # concatenate them into a Pipe object
  search = gt.Pipe([covol, covar])


