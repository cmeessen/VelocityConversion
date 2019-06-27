# -*- coding: utf-8 -*-
from VelocityConversion import MantleConversion

MC = MantleConversion()

# Data
MC.FileIn = "./VsSL2013.dat"
MC.LoadFile()

"""
The data can also be directly provided as a numpy array of the shape (nrows, 4)
where the columns are [x, y, z, v]. In this case, the array has to be assigned
to DataRaw, for example:

>> SomeArray = np.zeros([10,4])
>> MC.DataRaw = SomeArray
"""

# Set Velocity type
MC.SetVelType("S")

"""
Mineralogy
==========

The mineralogy can be read from a file, e.g.

>> MC.LoadMineralogy("ExampleMineralogy.csv")

or be set manually by specifying a dictionary with phases and fractions
"""
assemblage = {
    "ol": 0.617,
    "cpx": 0.133,
    "opx": 0.052,
    "gnt": 0.153,
    "jd": 0.045,
    "XFe": 0.11
}
MC.SetMineralogy(assemblage)

"""
Optional parameters
===================

Iron content
------------
The iron content can also be defined manually afterwars using

>> MC.SetXFe(0.1)

Attenuation mode
----------------
The attenuation mode can be assigned using

>> MC.SetQMode()

Thermal expansion coefficient
-----------------------------
The P/T dependency can be assigned using

>> MC.SetAlpha(var)

where var is a str with `const`, `T` or `PT`.
"""

# Create the tables containing velocities, temperatures and densities
MC.FillTables()

# Run the interpolation
MC.CalcPT()

# The results are stored in arrays named `Result_T` and `Result_Rho`. Their
# order corresponds to the same order as in `DataRaw`
MC.Result_T
MC.Result_Rho

# By specifying a filename, the tables can also be saved
MC.SaveFile('VsSL2013_module.dat')

