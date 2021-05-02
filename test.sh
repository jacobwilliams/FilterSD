#!/bin/bash

#
# Just a simple test of the individual routines.
#

rm -rf bin
mkdir bin

# COMPILER=gfortran
# FLAGS="-O2 -std=legacy"

COMPILER=ifort
FLAGS="-O0 -std18"

# To solve an NLP using a sparse matrix data structure, the subroutines in the following files are required
$COMPILER $FLAGS ./src/filterSD.f90 ./src/checkd.f90 ./src/glcpd.f90 ./src/l1sold.f90 ./src/shared.f90 ./src/schurQR.f90 ./src/sparseA.f90 ./src/util.f90 ./src/tests/hs106.f90 -o ./bin/hs106

# To solve an NLP using a dense matrix data structure, the subroutines in the following files are required
$COMPILER $FLAGS ./src/filterSD.f90 ./src/checkd.f90 ./src/glcpd.f90 ./src/l1sold.f90 ./src/shared.f90 ./src/denseL.f90 ./src/denseA.f90 ./src/util.f90 ./src/tests/hs106d.f90 -o ./bin/hs106d

# To solve an LCP using a sparse matrix data structure, the subroutines in the following files are required
$COMPILER $FLAGS ./src/glcpd.f90 ./src/shared.f90 ./src/util.f90 ./src/checkg.f90 ./src/schurQR.f90 ./src/sparseA.f90 ./src/tests/hs72.f90 -o ./bin/hs72

# To solve an LCP using a dense matrix data structure, the subroutines in the following files are required
$COMPILER $FLAGS ./src/glcpd.f90 ./src/checkg.f90 ./src/shared.f90 ./src/denseL.f90 ./src/denseA.f90 ./src/util.f90 ./src/tests/hs72d.f90 -o ./bin/hs72d
