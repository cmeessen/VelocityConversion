Minimum working example
=======================

This is the minimum working example from the README::

    from VelocityConversion import MantleConversion
    MC = MantleConversion()
    MC.LoadFile("./Examples/VsSL2013.dat")
    MC.SetVelType("S")
    MC.DefaultMineralogy()
    MC.FillTables()
    MC.CalcPT()
    MC.SaveFile("./Examples/VsSL2013_out.dat")

Functions used in the working example
-------------------------------------

.. autofunction:: VelocityConversion.MantleConversion.LoadFile

.. autofunction:: VelocityConversion.MantleConversion.SetVelType

.. autofunction:: VelocityConversion.MantleConversion.DefaultMineralogy

.. autofunction:: VelocityConversion.MantleConversion.FillTables

.. autofunction:: VelocityConversion.MantleConversion.CalcPT

.. autofunction:: VelocityConversion.MantleConversion.SaveFile

Optional settings
-----------------

.. autofunction:: VelocityConversion.MantleConversion.SetMineralogy

.. autofunction:: VelocityConversion.MantleConversion.SetXFe

.. autofunction:: VelocityConversion.MantleConversion.SetQMode

.. autofunction:: VelocityConversion.MantleConversion.SetAlpha
