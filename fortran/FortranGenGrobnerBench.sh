#!/bin/bash
set -e

cd "$(dirname "$0")"

echo "Compiling Fortran Generic Grobner Tests..."
gfortran -O3 -cpp -I./generic -J./generic \
    generic/IntModP.f90 \
    generic/VecExponent.f90 \
    generic/BitPackedExponent.f90 \
    generic/GenGrobner.f90 \
    generic/TestGenGrobner.f90 \
    -o generic/TestGenGrobner

echo "Compilation successful. Running benchmarks..."

echo "--------------------------------"
echo "Cyclic 4 (Vec Exponent)"
time ./generic/TestGenGrobner 4 0

echo "--------------------------------"
echo "Cyclic 4 (Bit-Packed Exponent)"
time ./generic/TestGenGrobner 4 1

echo "--------------------------------"
echo "Cyclic 5 (Vec Exponent)"
time ./generic/TestGenGrobner 5 0

echo "--------------------------------"
echo "Cyclic 5 (Bit-Packed Exponent)"
time ./generic/TestGenGrobner 5 1

echo "--------------------------------"
echo "Done."
