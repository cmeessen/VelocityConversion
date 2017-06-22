# VelocityConversion

This code is a python implementation of the p- and s-wave velocity to density
conversion approach after Goes et al. (2000). The implementation was optimised
for regular 3D grids using lookup tables instead of Newton iterations.

Goes et al. (2000) regard the expansion coefficient as temperature dependent
using the relation by Saxena and Shen (1992). In `Conversion.py`, the user can
additionally choose between a constant expansion coefficient or a pressure- and
temperature dependent coefficient that was derived from Hacker and Abers (2004).

For detailed information on the physics behind the approach have a look at the
original paper by Goes et al. (2000).

## Recommended citation for VelocityConversion
Meeßen, Christian (2017): VelocityConversion. V. v1.0.1. GFZ Data Services.
http://doi.org/10.5880/GFZ.6.1.2017.001

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

## How to get started

To get help on the usage type

```
python Conversion.py --help
```

The steps to prepare a conversion are

- definition of mantle rock composition in a `*.csv` file using the mineral
  terminology of `MinDB.csv`
- provide a velocity distribution on a regular 3D grid where columns are `x y z
  v`
- run the `Conversion.py` specifying the velocity type with `-type P` or
  `-type S`

Working examples for conversions are given in the `./Example/` directory.

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
