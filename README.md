# VelocityConversion

Conversion of P and S wave velocities to temperature and density after the
approach by Goes et al. (2000). The method is optimised to work with 3D regular
grids.

## How it works in general:

- define mantle rock composition
- load velocity distribution
- analyse number of depth values in the velocity distribution
- calculate synthetic temperature and density for the given composition and
  depths in a predefined temperature range, e.g. 300 to 3000 K
- for every point in the velocity distribution compare the velocity with the
  synthetic values and obtain density and temperature

To get started open a console and type

```
python Conversion.py --help
```

Examples for conversions are given in the `./Example/` directory.

## References

Berckhemer, H., W. Kampfmann, E. Aulbach, and H. Schmeling. “Shear Modulus and
Q of Forsterite and Dunite near Partial Melting from Forced-Oscillation
Experiments.” Physics of the Earth and Planetary Interiors, Special Issue
Properties of Materials at High Pressures and High Temperatures, 29, no. 1
(July 1, 1982): 30–41. doi:10.1016/0031-9201(82)90135-2.

Goes, S., R. Govers, and P. Vacher. “Shallow Mantle Temperatures under Europe
from P and S Wave Tomography.” Journal of Geophysical Research 105, no. 11
(2000): 153–11.

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
