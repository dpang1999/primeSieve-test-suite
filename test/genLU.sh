TIMEFORMAT='%3R'
results_file="results/genLU_results.txt"
sanity_file="results/sanity$1.txt"
#> "$results_file"
> "$sanity_file"
#Echo parameter given to script
echo "Running tests for genLU up to $1 x $1..." | tee -a "$sanity_file"

# Run Java test
echo -n "Java,$1," | tee -a "$results_file" "$sanity_file"
cd ../java/generic && javac *.java > /dev/null 2>&1
{ time java -cp ./ genLU $1 >> ../test/"$sanity_file" 2>&1;} 2>> ../test/"$results_file"
cd ../test