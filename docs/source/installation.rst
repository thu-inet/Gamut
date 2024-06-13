Installation
==================


GAMUT can be easily installed using **pip** with the following shell command:

.. code-block:: bash

    pip install gamut

Or you can clone GAMUT from its `source code <https://github.com/thu-inet/Gamut>`_ with the following shell command:

.. code-block:: bash

    git clone https://github.com/thu-inet/Gamut.git


The most recommended way is to create an isolated Python enrivoment for GAMUT using **Conda**. This can be done with the following shell command:

.. code-block:: bash

    conda create -n gamut python=3.10
    conda activate gamut
    conda install -c conda-forge gamut

After installation, you can use GAMUT in **Jupyter Notebook** with the following command:

.. code-block:: bash

    conda activate gamut
    jupyter notebook
