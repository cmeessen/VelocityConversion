#!/bin/bash

# Using constant expansion coefficient
python ../Conversion.py VsSL2013.dat -comp ExampleMineralogy.csv -type S \
-out VsSL2013_AlphaConst.dat

# Using T-dependent expansion coefficient
python ../Conversion.py VsSL2013.dat -comp ExampleMineralogy.csv -type S \
-AlphaT -out VsSL2013_AlphaT.dat

# Using P-T-dependent expansion coefficient
python ../Conversion.py VsSL2013.dat -comp ExampleMineralogy.csv -type S \
-AlphaPT -out VsSL2013_AlphaPT.dat
