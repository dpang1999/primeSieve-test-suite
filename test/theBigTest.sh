# Number of times to run the script
X=$1
#echo "$X"

warm_up=1000
./primeSieve.sh "$warm_up"

# Clear the results directory
rm -rf /home/vscode/primeSieve-test-suite/test/results/*


# Loop to run primeSieve.sh X times with different parameters
for ((i=2; i<=X+1; i++))
do
    param_value=$(eval echo \$$i)
    #echo "$param_value"
    ./primeSieve.sh "$param_value"
done