echo "FFT tests"
cd ../java && javac $(find . -name "*.java") > /dev/null 2>&1
echo -n "Java Specialized,1048576,"
time java -Xmx27G specialized.FFT 1048576
echo -n "Java Specialized,16777216,"
time java -Xmx27G specialized.FFT 16777216
echo -n "Java Specialized,67108864,"
time java -Xmx27G specialized.FFT 67108864
echo -n "Java Specialized,268435456,"
time java -Xmx27G specialized.FFT 268435456
cd ../test