# Run Go test
echo "Running Go test..."
go run ../src/primeSieve.go 100 > /dev/null 2>&1

# Run Julia test
echo "Running Julia test..."
julia ../src/primeSieve.jl 100 > /dev/null 2>&1

# Run Java test
echo "Running Java test..."
javac ../src/PrimeSieve.java > /dev/null 2>&1
java -cp ../src PrimeSieve 100 > /dev/null 2>&1

# Run TypeScript test
echo "Running TypeScript test..."
tsx ../src/index.ts 100 > /dev/null 2>&1

# Run C++ test
echo "Running C++ test..."
g++ -o ../src/primeSieve ../src/primeSieve.cpp > /dev/null 2>&1
../src/primeSieve 100 > /dev/null 2>&1

# Run Rust test
echo "Running Rust test..."
(cd ../rust && cargo run -- 100 > /dev/null 2>&1)

# Run C# test
echo "Running C# test..."
(cd ../csharp && dotnet run -- 100 > /dev/null 2>&1)

# Run Fortran test
echo "Running Fortran test..."
gfortran ../src/primeSieve.f90 -o ../src/primeSieve > /dev/null 2>&1
../src/primeSieve 100 > /dev/null 2>&1