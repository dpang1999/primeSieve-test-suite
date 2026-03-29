#!/bin/bash

# Navigate to the directory containing the script
cd "$(dirname "$0")"

echo "Compiling C++ Grobner implementations..."
g++ -O3 -std=c++17 generic/GenGrobner.cpp -o generic/GenGrobner
g++ -O3 -std=c++17 specialized/FiniteGrobner.cpp -o specialized/FiniteGrobner
g++ -O3 -std=c++17 specialized/GrobnerSmart.cpp -o specialized/GrobnerSmart
echo "Compilation finished."
echo "----------------------------------------"
time ./specialized/FiniteGrobner 4
time ./generic/GenGrobner 4 0
time ./specialized/FiniteGrobner 5
time ./generic/GenGrobner 5 0
time ./specialized/FiniteGrobner 6
time ./generic/GenGrobner 6 0

time ./specialized/GrobnerSmart 4
time ./generic/GenGrobner 4 1
time ./specialized/GrobnerSmart 5
time ./generic/GenGrobner 5 1
time ./specialized/GrobnerSmart 6
time ./generic/GenGrobner 6 1
