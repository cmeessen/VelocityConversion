.. _as_module:

Python module
=============

Follow the steps outlined in `_installation`, in order to be able to use
:mod:`VelocityConversion` as a Python module. Start with importing the main
class:

.. ipython:: python

    from VelocityConversion import MantleConversion
    MC = MantleConversion()

Input data
----------

The input data can be loaded from a **file** or provided as a **numpy array**.
The data should be organised in columns:

    +------------+----------+-----------------+----------+
    | **Column** | **Name** | **Description** | **Unit** |
    +------------+----------+-----------------+----------+
    | 0          | X        | A coordinate    | Anything |
    +------------+----------+-----------------+----------+
    | 1          | Y        | A coordinate    | Anything |
    +------------+----------+-----------------+----------+
    | 2          | Z        | Depth           | masl     |
    +------------+----------+-----------------+----------+
    | 3          | v        | Velocity        | m/s      |
    +------------+----------+-----------------+----------+

To load data from a file, use
:meth:`~VelocityConversion.MantleConversion.LoadFile`.

.. ipython:: python

    MC.LoadFile('../../Examples/VsSL2013.dat')


Alternatively, the data can be provided as a numpy array using
:meth:`~VelocityConversion.MantleConversion.LoadArray`. When providing the data
as a numpy array, it should have the shape ``[nrows, 4]``. This example shows
the numpy array structure by simply displaying the array loaded from the file
above:

.. ipython:: python

    # Adjust numpy just for some prettier printing of the array
    np.set_printoptions(precision=2, suppress=True)
    MC.DataRaw

In any case, the velocity type needs to be provided with
:meth:`~VelocityConversion.MantleConversion.SetVelType` (``S`` or ``P``):

.. ipython:: python

    MC.SetVelType('S')

Defining a mantle assemblage
----------------------------

The mineralogical assemblage of the mantle rocks must be provided. The
assemblage is defined by denoting the fractional proportion of a mineral phase
of the mantle rock. It can be defined either by using a Python dictionary, or
loaded from a file using :meth:`~VelocityConversion.LoadMineralogy`:

.. ipython:: python

    assemblage = {
        "ol": 0.67,
        "cpx": 0.045,
        "opx": 0.225,
        "gnt": 0.06,
        "XFe": 0.11
    }

.. note::

    ``XFe`` is the molar iron content of the rock and is optional. It can also
    be defined by calling :meth:`~VelocityConversion.MantleConversion.SetXFe`.

The available mineral phases are provided in the ``MinDB.csv``. This file can
be edited freely if new phase information is available.

To defined the composition, call
:meth:`~VelocityConversion.MantleConversion.SetMineralogy`:

.. ipython:: python

    MC.SetMineralogy(assemblage)

or call :meth:`~VelocityConversion.MantleConversion.DefaultMineralogy` in order
to set the default composition (garnet lherzolite by Jordan, 1979).


Start the conversion
--------------------

To start the conversion, simply run
:meth:`~VelocityConversion.MantleConversion.Convert`:

.. ipython:: python

    MC.Convert()

Results
-------

The results are stored in two 1D arrays named ``Result_T`` and ``Result_Rho``.
The order of the data in these arrays corresponds to the order of the points
in the input file.

.. ipython:: python

    print(MC.Result_T)
    print(MC.Result_Rho)

.. warning::
    The temperatures in ``Result_T`` are in **Kelvin**!

We can plot the results for 50km depth:

.. ipython::

    In [1]: import matplotlib.pyplot as plt
       ...: from matplotlib.tri import Triangulation
       ...:
       ...: def make_plot(MCObject, depth):
       ...:     fig, ax = plt.subplots(1, 3, sharex=True, sharey=True)
       ...:     indices = MCObject.DataRaw[:, 2] == depth
       ...:     depth_slice = MCObject.DataRaw[indices]
       ...:     velocity = depth_slice[:, 3]
       ...:     temperature = MCObject.Result_T[indices]
       ...:     density = MCObject.Result_Rho[indices]
       ...:     tri = Triangulation(depth_slice[:, 0]/1000.,
       ...:                         depth_slice[:, 1]/1000.)
       ...:     m_v = ax[0].tricontourf(tri, velocity, levels=np.arange(3500, 4700, 100))
       ...:     m_t = ax[1].tricontourf(tri, temperature, cmap='hot', levels=np.arange(950, 1650, 50))
       ...:     m_d = ax[2].tricontourf(tri, density, cmap='plasma', levels=np.arange(3280, 3340, 5))
       ...:     for a in ax:
       ...:         a.set_aspect('equal')
       ...:     fig.colorbar(m_v, ax=ax[0], label='Vs / m/s', orientation='horizontal',
       ...:                  ticks=[3500, 4000, 4500])
       ...:     fig.colorbar(m_t, ax=ax[1], label='T / K', orientation='horizontal',
       ...:                  ticks=[1100, 1300, 1500])
       ...:     fig.colorbar(m_d, ax=ax[2], label='Density / kg/m3', orientation='horizontal',
       ...:                  ticks=[3280, 3300, 3320, 3340])

    @savefig figure.png
    In [2]: make_plot(MC, -50e3)


Saving to file
~~~~~~~~~~~~~~~

The results can be saved to an output file with
:meth:`~VelocityConversion.MantleConversion.SaveFile`:

.. ipython:: python

    MC.SaveFile('fileout.dat')

.. testsetup:: *

    import os
    os.remove("fileout.dat")

The output file will contain metadata about the conversion parameters, for
example:

.. code-block:: none

    # Temperature output
    # Input file: VsSL2013.dat
    # Velocity scale factor: 1.0
    # Mantle composition:
    # cpx - 0.133
    # gnt - 0.153
    # jd - 0.045
    # ol - 0.617
    # opx - 0.052
    # XFe - 0.11
    # Pressure calculation: AK135
    # Wave frequency (Omega) / Hz: 0.02
    # Anelasticity parameters: Sobolev et al. (1996)
    # Alpha depending on: Nothing
    # Columns:
    # 1 - X
    # 2 - Y
    # 3 - Z / masl
    # 4 - V_S / m/s
    # 5 - T_syn / degC
    # 6 - Rho / kg/m3

.. _optional_settings:

Optional settings
-----------------

.. _pressure_reference_model:

Pressure reference model
~~~~~~~~~~~~~~~~~~~~~~~~

Pressure in :mod:`VelocityConversion` is computed one-dimensionally using the
earth reference model `AK135
<http://rses.anu.edu.au/seismology/ak135/ak135f.html>`_. If desired, pressure
calculation can be performed using a homogeneous density. Therefore, define

.. ipython:: python

    MC.SimpleP = True
    MC.SimpleRho = 3000.  # kg/m3, the default value

If pressure was computed in this way, a note in the output metadata will be
added:

.. code-block:: none

    # Pressure calculation: Simplified
    # Pressure calculation density: 3000.0 kg/m3

.. _attenuation_parameters:

Attenuation parameters
~~~~~~~~~~~~~~~~~~~~~~~

The attenuation model used can be changed. By default, parameterisation after
Sobolev et al. (1996) is used. If desired, it can be changed to the parameters
by Berckhemer et al. (1982) using
:meth:`~VelocityConversion.MantleConversion.SetQMode`:

.. ipython:: python

    MC.SetQMode(1)  # Sobolev et al. (1996)
    MC.SetQMode(2)  # Berckhemer et al. (1982)

.. _thermal_expansion_coefficient:

Thermal expansion coefficients
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Goes et al. (2000) use a temperature-dependent isobaric thermal expansion
coefficient. In general, the expansion coefficient will increase with
increasing pressure. However, the expansion coefficient is also sensitive to
pressure changes such that it will decrease with increasing pressure.
Temperature and pressure act against each other. I therefore implemented three
different options to use this expansion coefficient. They are defined using
:meth:`~VelocityConversion.MantleConversion.SetAlpha`:

+----------+-----------------------------------------------------------------+
| ``mode`` | **Explanation**                                                 |
+----------+-----------------------------------------------------------------+
| const    | Constant. Use :math:`\alpha_0` of Saxena and Shen (1992)        |
+----------+-----------------------------------------------------------------+
| T        | T-dependent. Use the formulation of Saxena and Shen (1992)      |
+----------+-----------------------------------------------------------------+
| | PT     | | Pressure- and temperature dependent. Uses data extracted from |
| |        | | Hacker and Abers (2006)                                       |
+----------+-----------------------------------------------------------------+

The temperature dependent alpha is defined by

.. math::

    \alpha(T) = \alpha_0 + \alpha_1 T + \alpha_2 T^{-1} + \alpha_3 T^{-2}

Where the parameters :math:`\alpha_i` are taken from Saxena and Shen (1992).
They are provided in the file ``MinDB.csv``.

The data for pressure and temperature dependent behaviour were extracted from
the Excel worksheet by Hacker and Abers (2006) using a
`VBA script <https://github.com/cmeessen/VelocityConversion/blob/master/AlphaPT/AlphaPT.bas>`_.

The following figure shows how the expansion coefficients compare for
clinopyroxene (*cpx*):

.. image:: extra/alpha.png

.. ipython:: python

    MC.SetAlpha('const')

