.. _guide_define_a_pipenet:

Define a PipeNet
================


PipeNet is the superior class in GAMUT. It is used to represent an entire workflow for gamma
spectrum analysis tasks. It is the container for Node and Flow objects.


A Node object is a dictionary to store multiple Spectrum objects. Keys can be of any hashable
type such as the integer or string, and None is used as the placeholder.


.. code-block:: python

  # define a node containing three spectra
  spec1 = gt.Spectrum(np.random.rand(100))
  spec2 = gt.Spectrum(np.random.rand(100))
  # None servers as a placeholder
  node = gt.Node({1: spec1, '2':spec2, 5: None})

A Flow object is an Operator object decorated with the input/ouput information. Notice that
when an Operator object takes multiple of spectra as its input, the thrid argument should be a
list of length that eqauls to its inp_num.


.. code-block:: python

  # define two Operator objects
  wavelet = gt.TranslationInvarianceWaveletSmoother('dmey', 'quadratic-soft', order=2)
  strp = gt.OtherOperator.SNIPStripper(high_order=False)

  # input the Node 'input' and the specrtum of key==0
  # output into the Node 1 and the spectrum of key=='smoothed'
  fsmooth = gt.Flow(wavelet, 'input', 0, 1, 'smoothed') # 1st Node

  # input the Node 1 and the specrtum of key=='smoothed'
  # output into the Node 'output' and the spectrum of key==999
  fstrip = gt.Flow(strp, 1, 'smoothed', 'output', 999)

  # when an Operator takes multiple spectra as the input
  # input the Node 0 and specrta of key=='measured' and 'background'
  # output into the Node 1 and the spectrum of key==999
  fstrip_baseline = gt.Flow(gt.Stripper(), 0, ['measured', 'background'], 1, 'stripped')


When defining a PipeNet object, only Operator and Flow objects should be explicitly defined.
And the Node objects will be automatically created in the PipeNet object with placeholders.
Besides, one must assign a entry Flow using the entry attribute. The default value is 0.

.. code-block:: python

  # define a PipeNet and its entry label
  pipenet = gt.PipeNet([fsmooth, fstrip])
  pipenet.entry = 'input'

When the PipeNet object is called, it will execute the entry Flow, and send its output spectrum
into the correct place. Then, it will check if there’s any other executable Flow objects (i.e. all
input spectra are well defined rather than placeholders). If so, they will be executed according to
the input/output order. The PipeNet object will repeat the check until no Flow objects can be
executed. Finally, a dictionary containing all internal Node objects is returned, and all analysis
results can be retrived from the dictionary.

.. code-block:: python

  # run the PipeNet and retrive the results
  nodes = pipenet(input)
  node = nodes['output'] # retrive a Node
  analyzed = node[0]

