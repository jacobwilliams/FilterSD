#!/bin/bash

#
# Just a simple test of the individual routines.
#

rm -rf bin
# rm -rf build
mkdir bin
# mkdir build

COMPILER=gfortran
FLAGS="-O2 -std=legacy"

# To solve an NLP using a sparse matrix data structure, the subroutines in the following files are required
$COMPILER $FLAGS ./src/filterSD.f90 ./src/checkd.f90 ./src/glcpd.f90 ./src/l1sold.f90 ./src/shared.f90 ./src/schurQR.f90 ./src/sparseA.f90 ./src/util.f90 ./src/hs106.f90 -o ./bin/hs106

$COMPILER $FLAGS ./src/filterSD.f90 ./src/checkd.f90 ./src/glcpd.f90 ./src/l1sold.f90 ./src/shared.f90 ./src/schurQR.f90 ./src/sparseA.f90 ./src/util.f90 ./src/hs72.f90 -o ./bin/hs72

# To solve an NLP using a dense matrix data structure, the subroutines in the following files are required
$COMPILER $FLAGS ./src/filterSD.f90 ./src/checkd.f90 ./src/glcpd.f90 ./src/l1sold.f90 ./src/shared.f90 ./src/denseL.f90 ./src/denseA.f90 ./src/util.f90 ./src/hs106d.f90 -o ./bin/hs106d

$COMPILER $FLAGS ./src/filterSD.f90 ./src/checkd.f90 ./src/glcpd.f90 ./src/l1sold.f90 ./src/shared.f90 ./src/denseL.f90 ./src/denseA.f90 ./src/util.f90 ./src/hs72d.f90 -o ./bin/hs72d


# # To solve an LCP using a sparse matrix data structure, the subroutines in the following files are required
# gfortran -O2 -std=legacy -o ./src/glcpd.f90 ./src/shared.f90 ./src/util.f90 ./src/checkg.f90 ./src/schurQR.f90 ./src/sparseA.f90  
# ar cr ./build/lcp_sparse.a *.o
# rm *.o

# # To solve an LCP using a dense matrix data structure, the subroutines in the following files are required
# gfortran -O2 -std=legacy -o ./src/glcpd.f90 ./src/checkg.f90 ./src/shared.f90 ./src/denseL.f90 ./src/denseA.f90 ./src/util.f90
# ar cr ./build/lcp_dense.a *.o
# rm *.o

# # test programs:
# gfortran -O2 -std=legacy ./src/hs106.f90   ./build/nlp_sparse.a -o ./bin/hs106
# gfortran -O2 -std=legacy ./src/hs106d.f90  ./build/nlp_dense.a  -o ./bin/hs106d
# gfortran -O2 -std=legacy ./src/hs72.f90    ./build/nlp_sparse.a -o ./bin/hs72
# gfortran -O2 -std=legacy ./src/hs72d.f90   ./build/nlp_dense.a  -o ./bin/hs72d
