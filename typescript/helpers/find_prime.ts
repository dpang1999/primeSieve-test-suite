// Returns the smallest prime p such that p ≡ 1 mod n
import { primeSieve } from '../helpers/primeSieve';
export function find_prime_congruent_one_mod_n(n: number): number {
	let k = 1;
	while (true) {
		const p = k * n + 1;
		const primes = primeSieve(p + 1);
		if (!primes[p]) return p;
		k++;
	}
}
