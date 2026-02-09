package helpers;

import java.util.Arrays;
import specialized.primeSieve;

public class FindPrime {

    // Find the smallest prime p such that p ≡ 1 mod n
    public static int findPrimeCongruentOneModN(int n) {
        int k = 1;
        while (true) {
            int p = k * n + 1;
            boolean[] primes = primeSieve.sieve(p+1);
            if (!primes[p]) {
                return p;
            }
            k++;
        }
    }

    // Example usage
    public static void main(String[] args) {
        int n = 7;
        int p = findPrimeCongruentOneModN(n);
        System.out.println("Smallest prime ≡ 1 mod " + n + ": " + p);
    }
}
