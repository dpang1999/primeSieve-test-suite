javac $(find . -name "*.java") 
time java -Xmx27G specialized.LU 3000;
time java -Xmx27G specialized.LU 4000;
time java -Xmx27G specialized.LU 5000;
time java -Xmx27G generic.GenLU 3000 0;
time java -Xmx27G generic.GenLU 4000 0;
time java -Xmx27G generic.GenLU 5000 0;
time java -Xmx27G specialized.FiniteLU 3000;
time java -Xmx27G specialized.FiniteLU 4000;
time java -Xmx27G specialized.FiniteLU 5000;
time java -Xmx27G generic.GenLU 3000 1;
time java -Xmx27G generic.GenLU 4000 1;
time java -Xmx27G generic.GenLU 5000 1;