Getting Started
==================

Interactive Analysis with Basic APIs
-------------------------------------------------

Here is a MWE (minimum working example) of GAMUT to analyze a gamma spectrum interactively.

.. code-block:: python

    import gamut as gt

    # import a gamma spectrum
    spec = gt.Spectrum.from_GammaVision("spectrum.spe")

    # initiate an operator/algorithm
    smoother = gt.CentroidSmoother(order=2)

    # analyze the spectrum
    output = smoother(spec)

    # visualize the spectrum
    ax = output.plot()
    plt.gcf().savefig("smoothed_spectrum.png")

    # export the spectrum
    out.export_to_GammaVision("smoothed_spectrum.png")

The above code snippet does the following:

.. code-block:: python

    import gamut as gt

0. Import the GAMUT library. This is the first step of using GAMUT. The library is constructed in a hierarchical manner, with the top-level module being :py:mod:`gamut`. For more information on the available submodules, see :ref:`API References <api_reference>`.

.. code-block:: python

    # import a gamma spectrum
    spec = gt.Spectrum.from_GammaVision("spectrum.spe")

1. Import a gamma spectrum from a file. Here a :py:class:`Spectrum` object is initialized by importing from a GammaVision file. For more information on how to import a spectrum from different file formats, see :ref:`Import a Spectrum <guide_import_a_spectrum>`.

.. code-block:: python

    # import a gamma spectrum
    smoother = gt.CentroidSmoother(order=2)

2. Initialize an operator/algorithm. Here, a :py:class:`CentroidSmoother` object is initialized with a smoothing order of 2. Parameters are passed to the algorithm during initialization. For more information on the available algorithms and their initialization, see :ref:`List of Algorithms <algotithms>` and :ref:`Initialize an Operator <guide_initialize_an_operator>`.


.. code-block:: python

    # analyze the spectrum
    output = smoother(spec)

3. Analyze the spectrum. One must pass the correct number of :py:class:`Spectrum` objects to the algorithm (specified by :py:attr:`inp_num` attribute and it is 1 for :py:class:`CentroidSmoother`). And a :py:class:`Spectrum` object is returned as output.

.. code-block:: python

    # visualize the spectrum
    ax = output.plot()
    spec.plot(ax)
    plt.gcf().savefig("smoothed_spectrum.png")

4. Visualize the spectrum. Multiple visualization methods are available to control which data to display. For more information on how to visualize a spectrum, see :ref:`Visualize a Spectrum <guide_visualize_a_spectrum>`.

.. code-block:: python

    # export the spectrum
    out.export_to_GammaVision("smoothed_spectrum.png")

5. Export the spectrum. The **Spectrum** object can be exported to a file. For more information on how to export a spectrum to different file formats, see :ref:`Export a Spectrum <guide_export_a_spectrum>`.

End-to-End Analysis with the Predefined Workflow
------------------------------------------

There is also a pre-defined workflow in GAMUT, which provides the out-of-the-box gamma spectrum analysis capacity for beginner users.

.. code-block:: python

    import gamut as gt

    # import a gamma spectrum
    spec = gt.Spectrum.from_GammaVision("spectrum.spe")

    # run the default workflow
    nodes = gt.default_workflow(spec)

    # analyze the spectrum
    output = nodes[4][0]

    # export the spectrum
    out.export_to_GammaVision("analyzed_spectrum.png")

The key part of the above code snippet does the following:

.. code-block:: python

    # run the default workflow
    nodes = gt.default_workflow(spec)

1. Run the default workflow. GAMUT provides a default workflow for the out-of-the-box gamma spectrum analysis. It is a :py:class:`PipeNet` object, and can execute a serie of algorithms automatically upon calling. For more information on the :py:class:`PipeNet` class, see :ref:`PipeNet <guide_define_a_pipenet>`.

.. code-block:: python

    # run the default workflow
    output = nodes[4][0]

1. Retrive the analyzed spectrum. After executing the :py:class:`PipeNet` object, the output is a list of :py:class:`Node` objects, which is used to store :py:class:`Spectrum` objects. In the default workflow, the final analyzed spectrum is stored in the 5th node. So we retrieve the first element in the fifth :py:class:`Node` object to get the analysis results.


