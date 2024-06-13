.. _guide_initialize_an_operator:

Initialize an Operator
======================

:py:class:`Operator` is another key class in GAMUT. All algorithms are defined as subclasses of Operator in
GAMUT, and they must be initialized with configuration parameters before applying to gamma spectra.

.. code-block:: python

  # initialize a Savitzky-Golay smoothing algorithm
  savit = gt.SavitzkySmoother(order=2, hwidth=3)
  savit.order = 3 # parameters can be modified later


Then :py:class:`Operator` objects can be called to apply the algorithms. It should be noted that the number
of input spectra must equals to its :py:attr:`inp_num` attribute, except when :py:attr:`inp_num` euqals to -1, which
means any number of input spectra. Most of :py:class:`Operator` objects’ inp_num objects equals to 1.
The number of output spectra is always 1.

.. code-block:: python

  # apply the Savitzky-Golay smoothing algorithm
  smoothed = savit(spectra)


If the :py:attr:`inp_num` of an :py:class:`Operator` object equals to 1, then the pipeline operator **»** can be used to
simplify the initilization and calling. It is also used to call multiple :py:class:`Operator` objects consecutively.

.. code-block:: python

  # input the spectrum into the operator
  smoothed = spec >> gt.SavitzkySmoother(order=2, hwidth=3)# assume that op1, op2, and op3 are three operators of inp_num=1
  analyzed = spec >> op1 >> op2 >> op3

:py:class:`Functionor` is a special subclass of :py:class:`Operator`, which can wrap any function into an Operator
object.

.. code-block:: python

  # define an operator that combines two spectra
  def combine(spec1, spec2):
  spec2.regions = deepcopy(spec1.regions)
  return spec2
  comb = gt.Functionor(inp_num=2, func=combine)


