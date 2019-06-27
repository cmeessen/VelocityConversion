#!/bin/bash
set -e
export PATH=$PWD:$PATH

# Using constant expansion coefficient
VelocityConversion VsSL2013.dat -comp ExampleMineralogy.csv -type S \
-out VsSL2013_AlphaConst.dat

# Using T-dependent expansion coefficient
VelocityConversion VsSL2013.dat -comp ExampleMineralogy.csv -type S \
-AlphaT -out VsSL2013_AlphaT.dat

# Using P-T-dependent expansion coefficient
VelocityConversion VsSL2013.dat -comp ExampleMineralogy.csv -type S \
-AlphaPT -out VsSL2013_AlphaPT.dat
