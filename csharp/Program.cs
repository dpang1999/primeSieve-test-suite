// See https://aka.ms/new-console-template for more information
class Program
{
    public static Boolean[] primeSieve(int num) {
        Boolean[] primes = new Boolean[num];
        // default instantiation is false, flipping semantics
        primes[0] = true;
        primes[1] = true;
        for (int i = 2; i < num; i++) {
            if (!primes[i]) {
                int j = i;
                while (j * i < num) {
                    primes[j*i] = true;
                    j++;
                }
            }
        }
        return primes;
    }

    static void Main(string[] args)
    {
        int max = 42;
        if (args.Length > 0) {
            max = int.Parse(args[0]);
        }
        Boolean[] temp = primeSieve(max);
        for (int i = 2; i<temp.Length; i++) {
            if (!temp[i]) {
                Console.WriteLine(i);
            }
        }
    }
}
