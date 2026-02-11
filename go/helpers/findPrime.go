package helpers

// FindCongruentPrime finds a prime number p such that p ≡ 1 (mod size)
// This is useful for FFT operations in finite fields where the prime
// must satisfy p = k*size + 1 for some integer k
func FindCongruentPrime(size int) int {
	if size == 0 {
		return 0
	}

	// Start searching from the first candidate
	candidate := size + 1

	for {
		if isPrime(candidate) {
			return candidate
		}
		candidate += size
	}
}

// isPrime checks if a number is prime using trial division
func isPrime(n int) bool {
	if n < 2 {
		return false
	}
	if n == 2 {
		return true
	}
	if n%2 == 0 {
		return false
	}

	// Trial division up to sqrt(n)
	i := 3
	for i*i <= n {
		if n%i == 0 {
			return false
		}
		i += 2
	}
	return true
}
