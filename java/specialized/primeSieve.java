package specialized;
public class primeSieve {
    public static boolean[] sieve(int num) {
        boolean[] primes = new boolean[num];
        // default instantiation is false, flipping semantics
        primes[0] = true;
        primes[1] = true;
        for (int i = 2; i < num; i++) {
            if (!primes[i]) {
                int j = i;
                while (j * i < num && j * i > 0) { // check for overflow
                    primes[j*i] = true;
                    j++;
                }
            }
        }
        return primes;
    }
    
    public static void main(String[] args) throws Exception {
        int max = 42;
        if (args.length > 0){
            max = Integer.parseInt(args[0]);
        }
        boolean[] temp = sieve(max);
        for (int i = 2; i < max; i++) {
            if (!temp[i]) {
                System.out.println(i);
            }
        }
    }
}
