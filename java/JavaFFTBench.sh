javac $(find . -name "*.java") 
time java -Xmx27G generic.GenFFT 1048576 0;
time java -Xmx27G generic.GenFFT 1048576 1;
time java -Xmx27G generic.GenFFT 16777216 0;
time java -Xmx27G generic.GenFFT 16777216 1;
time java -Xmx27G generic.GenFFT 67108864 0;
time java -Xmx27G generic.GenFFT 67108864 1;
time java -Xmx27G specialized.FFT 1048576;
time java -Xmx27G specialized.FFT 16777216;
time java -Xmx27G specialized.FFT 67108864;
time java -Xmx27G specialized.FiniteFFT 1048576;
time java -Xmx27G specialized.FiniteFFT 16777216;
time java -Xmx27G specialized.FiniteFFT 67108864;