echo "Gen Java FFT tests"
cd ../java && javac $(find . -name "*.java") > /dev/null 2>&1
echo -n "Java Generic,67108864,"
time java -Xmx27G generic.GenFFT 67108864 1
echo -n "Java Generic,1048576,"
time java -Xmx27G generic.GenFFT 1048576 1
echo -n "Java Generic,16777216,"
time java -Xmx27G generic.GenFFT 16777216 1
cd ../test