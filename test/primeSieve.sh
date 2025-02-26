> primeSieve_results.txt
> sanity.txt
#Echo parameter given to script
echo "Running tests for primes up to $1..." | tee -a primeSieve_results.txt sanity.txt

# I'm conflicted on whether to include the compilation step in the timing or not
# C# for example includes compilation and running in the dotnet command whereas in java has separate commands
# But my question becomes, do people care about compilation time?

# Run Go test
echo "Running Go test..." | tee -a primeSieve_results.txt sanity.txt
{ time go run ../go/primeSieve.go $1 >> ../test/sanity.txt 2>&1;} 2>> primeSieve_results.txt

# Run Julia test
echo "Running Julia test..." | tee -a primeSieve_results.txt sanity.txt
{ time julia ../julia/primeSieve.jl $1 >> ../test/sanity.txt 2>&1;} 2>> primeSieve_results.txt

# Run Java test
echo "Running Java test..." | tee -a primeSieve_results.txt sanity.txt
cd ../java && javac primeSieve.java > /dev/null 2>&1
{ time java -cp ./ primeSieve $1 >> ../test/sanity.txt 2>&1;} 2>> ../test/primeSieve_results.txt
cd ../test

# Run TypeScript test
echo "Running TypeScript test..." | tee -a primeSieve_results.txt sanity.txt
{ time tsx ../typescript/index.ts $1 >> ../test/sanity.txt 2>&1;} 2>> primeSieve_results.txt

# Run C++ test
echo "Running C++ test..." | tee -a primeSieve_results.txt sanity.txt
cd ../cpp && g++ -o primeSieve primeSieve.cpp > /dev/null 2>&1
{ time ./primeSieve $1 >> ../test/sanity.txt 2>&1;} 2>> ../test/primeSieve_results.txt
cd ../test

# Run Rust test
echo "Running Rust test..." | tee -a primeSieve_results.txt sanity.txt
cd ../rust
{ time cargo run -- $1 >> ../test/sanity.txt 2>&1;} 2>> ../test/primeSieve_results.txt
cd ../test

# Run C# test
echo "Running C# test..." | tee -a primeSieve_results.txt sanity.txt
cd ../csharp
{ time dotnet run -- $1 >> ../test/sanity.txt 2>&1;} 2>> ../test/primeSieve_results.txt
cd ../test

# Run Fortran test
echo "Running Fortran test..." | tee -a primeSieve_results.txt sanity.txt
cd ../fortran && gfortran primeSieve.f90 -o primeSieve > /dev/null 2>&1
{ time ./primeSieve $1 >> ../test/sanity.txt 2>&1;} 2>> ../test/primeSieve_results.txt
cd ../test
