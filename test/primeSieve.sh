results_file="results/primeSieve_results$1.txt"
sanity_file="results/sanity$1.txt"

> "$results_file"
> "$sanity_file"
#Echo parameter given to script
echo "Running tests for primes up to $1..." | tee -a "$results_file" "$sanity_file"

# Run Go test
echo "Running Go test..." | tee -a "$results_file" "$sanity_file"
{ time go run ../go/primeSieve.go $1 >> ../test/"$sanity_file" 2>&1;} 2>> "$results_file"

# Run Julia test
echo "Running Julia test..." | tee -a "$results_file" "$sanity_file"
{ time julia ../julia/primeSieve.jl $1 >> ../test/"$sanity_file" 2>&1;} 2>> "$results_file"

# Run Java test
echo "Running Java test..." | tee -a "$results_file" "$sanity_file"
cd ../java && javac primeSieve.java > /dev/null 2>&1
{ time java -cp ./ primeSieve $1 >> ../test/"$sanity_file" 2>&1;} 2>> ../test/"$results_file"
cd ../test

# Run TypeScript test
echo "Running TypeScript test..." | tee -a "$results_file" "$sanity_file"
{ time tsx ../typescript/index.ts $1 >> ../test/"$sanity_file" 2>&1;} 2>> "$results_file"

# Run C++ test
echo "Running C++ test..." | tee -a "$results_file" "$sanity_file"
cd ../cpp && g++ -o primeSieve primeSieve.cpp > /dev/null 2>&1
{ time ./primeSieve $1 >> ../test/"$sanity_file" 2>&1;} 2>> ../test/"$results_file"
cd ../test

# Run Rust test
echo "Running Rust test..." | tee -a "$results_file" "$sanity_file"
cd ../rust
{ time cargo run -- $1 >> ../test/"$sanity_file" 2>&1;} 2>> ../test/"$results_file"
cd ../test

# Run C# test
echo "Running C# test..." | tee -a "$results_file" "$sanity_file"
cd ../csharp
{ time dotnet run -- $1 >> ../test/"$sanity_file" 2>&1;} 2>> ../test/"$results_file"
cd ../test

# Run Fortran test
echo "Running Fortran test..." | tee -a "$results_file" "$sanity_file"
cd ../fortran && gfortran primeSieve.f90 -o primeSieve > /dev/null 2>&1
{ time ./primeSieve $1 >> ../test/"$sanity_file" 2>&1;} 2>> ../test/"$results_file"
cd ../test
