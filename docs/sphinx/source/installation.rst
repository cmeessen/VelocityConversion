.. _installation:

Installation
============

:mod:`VelocityConversion` requires Python 3 and numpy. Install everything using
``pip`` by running

.. code-block:: bash

    pip install numpy velocityconversion

Using the latest version not on PyPI
------------------------------------

If you want to use the very latest version, or want to contribute, clone the
repository to your local hard drive:

.. code-block:: bash

    git clone https://github.com/cmeessen/VelocityConversion.git

change into the ``VelocityConversion`` directory, and install an editable
version:

.. code-block:: bash

    pip install -e .

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
