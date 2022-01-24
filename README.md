# VelocityConversion

[![DOI](https://zenodo.org/badge/87794116.svg)](https://zenodo.org/badge/latestdoi/87794116)

- [VelocityConversion](#velocityconversion)
  - [Introduction](#introduction)
  - [Getting started](#getting-started)
    - [Use the latest version not on PyPI](#use-the-latest-version-not-on-pypi)
  - [Usage as command line tool](#usage-as-command-line-tool)
  - [Usage as a Python module](#usage-as-a-python-module)
  - [Modifying physical properties of the minerals](#modifying-physical-properties-of-the-minerals)
  - [Contributing](#contributing)
  - [Citing](#citing)
  - [References](#references)
  - [Licence](#licence)

## Introduction

This code is a python implementation of the p- and s-wave velocity to density
conversion approach after Goes et al. (2000). The implementation was optimised
for regular 3D grids using lookup tables instead of Newton iterations.

Goes et al. (2000) regard the expansion coefficient as temperature dependent
using the relation by Saxena and Shen (1992). In `VelocityConversion`, the user
can additionally choose between a constant expansion coefficient or a pressure-
and temperature dependent coefficient that was derived from Hacker and Abers
(2004).

For detailed information on the physics behind the approach have a look at the
original paper by Goes et al. (2000).

## Getting started

`VelocityConversion` requires Python 3 and numpy. Install `numpy` and
`VelocityConversion` by running

```bash
pip install numpy velocityconversion
```

To uninstall `VelocityConversion`, run

```bash
pip uninstall velocityconversion
```

### Use the latest version not on PyPI

If you want to use the very latest version, or want to
[contribute](#contributing), clone the repository to you local hard drive:

```bash
git clone https://github.com/cmeessen/VelocityConversion.git
```

or, if you haven an [SSH key](https://docs.github.com/en/free-pro-team@latest/github/authenticating-to-github/generating-a-new-ssh-key-and-adding-it-to-the-ssh-agent)
associated to your account:

```bash
git clone git@github.com:cmeessen/VelocityConversion.git
```

To check whether everything is working run the tests

```bash
python test.py
```

If the output looks like this, everything is working fine:

```
test_vp_AlphaConst (__main__.TestVelocityConversion) ... ok
test_vs_AlphaConst (__main__.TestVelocityConversion) ... ok
test_vs_AlphaPT (__main__.TestVelocityConversion) ... ok
test_vs_AlphaT (__main__.TestVelocityConversion) ... ok

----------------------------------------------------------------------
Ran 4 tests in 1.633s

OK
```

## Usage as command line tool

In order to use the code as command line tool, add the `./Examples` directory
to your `PATH`, e.g. in your bash profile:

```bash
export PATH=/path/to/VelocityConversion/Examples:$PATH
```

Alternatively you can move the bash script
[VelocityConversion](./Examples/VelocityConversion) to a place that is within
your `PATH`. Now the bash script `VelocityConversion` can be executed:

```
VelocityConversion

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
        --version
```

The steps to prepare a conversion are

- definition of mantle rock composition in a `*.csv` file using the mineral
  terminology of [MinDB.csv](./VelocityConversion/MinDB.csv)
- provide a velocity distribution on a regular 3D grid where columns are `x y z
  v`
- run `VelocityConversion` specifying the velocity type with `-type P` or
  `-type S`

Working examples for the usage as command line tool are provided in the script
[RunExamples.sh](./Examples/RunExamples.sh).

## Usage as a Python module

VelocityConversion can also be imported as a Python module. Therefore, navigate
to the folder that contains your clone of the repository (and
[setup.py](./setup.py)) and execute

```bash
pip install -e .
```

Now, the module can be imported to Python:

```python
from VelocityConversion import MantleConversion
MC = MantleConversion()
```

A short working example for a conversion is:

```python
from VelocityConversion import MantleConversion
MC = MantleConversion()
MC.LoadFile("./Examples/VsSL2013.dat")
MC.SetVelType("S")
MC.DefaultMineralogy()
MC.FillTables()
MC.CalcPT()
MC.SaveFile("./Examples/VsSL2013_out.dat")
```

For a more complete documentation on how to use `VelocityConversion` as a Python
module please visit the
[documentation](https://cmeessen.github.io/VelocityConversion/).

## Modifying physical properties of the minerals

The database that contains the physical properties of the individual mineral
phases is stored in [MinDB.csv](./VelocityConversion/MinDB.csv).
Mineral parameters can be edited, or new minerals added. A new mineral phase
should then be referred to in the code or the assemblage file using the name
that was assigned in the `phase` column of `MinDB.csv`.

## Contributing

Please see [CONTRIBUTING.md](./CONTRIBUTING.md) if you want to contribute to
`VelocityConversion`.

## Citing

If you use this code, please consider citing it as

> Meeßen, Christian (2019): "VelocityConversion (v1.1.2)". Zenodo,
> http://doi.org/10.5281/zenodo.5897455.

or refer to [CITATION.cff](./CITATION.cff).

## References

Berckhemer, H., W. Kampfmann, E. Aulbach, and H. Schmeling. “Shear Modulus and
Q of Forsterite and Dunite near Partial Melting from Forced-Oscillation
Experiments.” Physics of the Earth and Planetary Interiors, Special Issue
Properties of Materials at High Pressures and High Temperatures, 29, no. 1
(July 1, 1982): 30–41. doi:10.1016/0031-9201(82)90135-2.

Goes, S., R. Govers, and P. Vacher. “Shallow Mantle Temperatures under Europe
from P and S Wave Tomography.” Journal of Geophysical Research 105, no. 11
(2000): 153–11. doi:10.1029/1999jb900300.

Hacker, Bradley R., and Geoffrey A. Abers. “Subduction Factory 3: An Excel
Worksheet and Macro for Calculating the Densities, Seismic Wave Speeds, and H2O
Contents of Minerals and Rocks at Pressure and Temperature.” Geochemistry,
Geophysics, Geosystems 5, no. 1 (January 1, 2004): Q01005.
doi:10.1029/2003GC000614.

Kennett, B. L. N., E. R. Engdahl, and R. Buland. “Constraints on Seismic
Velocities in the Earth from Traveltimes.” Geophysical Journal International
122, no. 1 (July 1, 1995): 108–24. doi:10.1111/j.1365-246X.1995.tb03540.x.

Saxena, Surendra K., and Guoyin Shen. “Assessed Data on Heat Capacity, Thermal
Expansion, and Compressibility for Some Oxides and Silicates.” Journal of
Geophysical Research: Solid Earth 97, no. B13 (Dezember 1992): 19813–25.
doi:10.1029/92JB01555.

Schaeffer, A. J., and S. Lebedev. “Global Shear Speed Structure of the Upper
Mantle and Transition Zone.” Geophysical Journal International 194, no. 1 (July
1, 2013): 417–49. doi:10.1093/gji/ggt095.

Sobolev, Stephan V., Hermann Zeyen, Gerald Stoll, Friederike Werling, Rainer
Altherr, and Karl Fuchs. “Upper Mantle Temperatures from Teleseismic Tomography
of French Massif Central Including Effects of Composition, Mineral Reactions,
Anharmonicity, Anelasticity and Partial Melt.” Earth and Planetary Science
Letters 139, no. 1–2 (März 1996): 147–63. doi:10.1016/0012-821X(95)00238-8.

## Licence

Licence: GNU General Public Licence, Version 3, 29 June 2007

Copyright (2017): Christian Meeßen, Potsdam, Germany

VelocityConversion is free software: you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by the Free
Software Foundation, either version 3 of the License, or (at your option) any
later version. VelocityConversion is distributed in the hope that it will be
useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public
License for more details. You should have received a cop y of the GNU General
Public License along with this program. If not, see
http://www.gnu.org/licenses/.
