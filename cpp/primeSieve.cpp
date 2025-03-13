#include <iostream>
#include <vector>
// g++ -o (output file name) (source file name)
// ./(output file name)
std::vector<bool> primeSieve(int num)
{
    std::vector<bool> primes(num, true);
    primes[0] = false;
    primes[1] = false;
    for(long long i = 2; i < num; i++)
    {
        if(primes[i] != false) {
            long long j = i;
            while (i*j < num) {
                primes[i*j] = false;
                j++;
            }
        }
    }
    return primes;
}

int main(int argc, char *argv[])
{
    int max = 42; // default value
    if (argc > 1) // program name is always the first argument
    {
        max = atoi(argv[1]);
    }

    std::vector<bool> temp = primeSieve(max);
    for(int i = 0; i < max; i++)
    {
        if(temp[i] != false)
        {
            std::cout << i << std::endl;
        }
    }
    return 0;
}