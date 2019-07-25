.. _installation:

Installation
============

Prerequisites
-------------

:mod:`VelocityConversion` requires Python 2.7 or 3, and a few modules.
I recommend to use Python 3 and `Anaconda <https://continuum.io>`_.

With Anaconda
~~~~~~~~~~~~~

Install the required modules with Anaconda:

.. code-block:: bash

    conda install -c conda-forge numpy

To be able to build the documentation from source, further packages are
required:

.. code-block:: bash

    conda install -c conda-forge ipython sphinx m2r sphinx_rtd_theme

With pip
~~~~~~~~~

or if you are using a Python version other than Anaconda:

.. code-block:: bash

    pip install numpy

and to be able to compile the documentation:

.. code-block:: bash

    pip install sphinx m2r sphinx_rtd_theme


Installing
----------

:mod:`VelocityConversion` can be used both as a command line tool and as a
Python module. In both cases, the first step is to create a clone or download
the repository.

.. code-block:: bash

    git clone https://github.com/cmeessen/VelocityConversion.git

and install it with

.. code-block::

    pip install .


Run tests
~~~~~~~~~~

To check whether everything is running correctly, enter the newly created
folder and run:

.. code-block:: bash

    python test.py

If the output looks like this, everything is working well:

.. code-block::

    test_vp_AlphaConst (__main__.TestVelocityConversion) ... ok
    test_vs_AlphaConst (__main__.TestVelocityConversion) ... ok
    test_vs_AlphaPT (__main__.TestVelocityConversion) ... ok
    test_vs_AlphaT (__main__.TestVelocityConversion) ... ok

    ----------------------------------------------------------------------
    Ran 4 tests in 1.633s

    OK

Usage as a Python module
------------------------

:mod:`VelocityConversion` can be used as a Python module. Therefore,
navigate to the folder that contains your clone of the repository (and
[setup.py](./setup.py)) and execute

.. code-block:: bash

    pip install -e .


Now, the module can be imported to Python:

.. code-block:: python

    from VelocityConversion import MantleConversion
    MC = MantleConversion()


Detailed instructions on the module usage a given in :ref:`as_module`.

Usage as a command line tool
----------------------------

:mod:`VelocityConversion` can also be run from the terminal. Please refer
to :ref:`terminal` for more information. Using it as a terminal application
does not require a pip installation.

