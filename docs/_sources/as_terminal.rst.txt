.. _terminal:

Terminal application
====================

:mod:`VelocityConversion` can also be used as a terminal application. Just
create a file named ``VelocityConversion`` (or any other name) and
add these lines:

.. code-block:: bash

    #!/bin/bash -e

    # Get the directory of this script
    DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

    python ../VelocityConversion/__init__.py $@

Make sure that the path points to the correct location. The file needs to be
executable:

.. code-block:: bash

    chmod +x VelocityConversion

.. note::

    If you want the script to be available everywhere, put it into a folder
    that is within your ``PATH``, or add it for example in you ``.bashrc``:
    ``export PATH=./path/to/VelocityConversion:$PATH``.

Executing the script will output

.. code-block:: none

    Usage: VelocityConversion FileIn -type <P|S> [optional args]
        Optional arguments:
            -AlphaT
            -AlphaPT
            -dT <val>
            -comp <Filename>
            -h | --help
            -NN
            -out <FileOut>
            -scaleV <value>
            -setQ <1|2>
            -v | -verbose
            -XFe <val>

Getting help
------------

Simply run ``VelocityConversion -h`` to get a more extensive help directly
within the terminal:

.. code-block:: none

    VelocityConversion --help

    Inverts seismic velocities for temperature and density based
    on input composition. By default the expansion coefficient is
    treated as pressure- and temperature-independent.

    Usage of VelocityConversion:

    Minimum requirements
    --------------------

    VelocityConversion FileIn -type <P|S> -comp <Filename>
        FileIn
            Input file name and path
        -type <P|S>
            Defines wave type P or S
        -comp <Filename>
            Comma-separated file containing mantle rock assemblage
        Example
            VelocityConversion Input.dat -type S -comp pyrolite.csv
            The output file will be InputOut.dat

    Input file requirements
    -----------------------

        Input file order and units must be: X Y Z V
        X Y - Coordinates in any units
        Z   - Depth in meters below or above sea level
        V   - Seismic velocity in m/s
        Can contain header lines if they are marked, i.e. with #

    Optional arguments
    ------------------

        -AlphaT
            Calculate Alpha based on Temperature after Saxena and
            Shen (1992). Default: const.

        -AlphaPT
            Calculate P/T-dependent Alpha based on excel worksheet
            from Hacker and Abers (2004). Default: const.

        -dT
            Changes the temperature increment in the P-T tables.
            Default increment is 1K.

        -NN
            Output file header information reduced to # of points

        -out <FileOut>
            Define the output path and/or file name
            Example: -out ../Output.dat

        -scaleV <value>
            Scale the velocity with the given value

        -setQ <1|2>
            Define attenuation parameters after 1 - Sobolev et al.
            (1996) or 2 - Berckhemer et al. (1982). Default: 1.
            Example: -setQ 2

        -v | -verbose
            Displays debugging messages.

        -XFe <XFe> | -xfe <XFe>
            Define iron content XFe in mole fractions. Default
            value XFe = 0.0
            Example: -XFe 0.1



Defining the mineralogy
-----------------------

The mineralogy must be provided as a comma-separated text file where the first
line denotes the columns ``phase`` and ``fraction``. The phase names correspond
to the ones in the ``MinDB.csv``:

.. code-block:: none

    phase,fraction
    ol,0.617
    cpx,0.133
    opx,0.052
    gnt,0.153
    jd,0.045
    XFe,0.11

Optional arguments
------------------

The optional arguments are well explained in the termial help. Please refer to
:ref:`optional_settings` for more information about the available pressure
computation mehtods, expansion coefficients and attenuation parameters.