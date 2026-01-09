use crate::helpers::prime_sieve::prime_sieve;

pub fn find_prime_congruent_one_mod_n(n: usize) -> usize {
    let mut k = 1;
    loop {
        let p = k * n + 1;
        let primes = prime_sieve(p+1);
        if primes[p] { // You can use your own is_prime or sieve
            return p;
        }
        k += 1;
    }
}