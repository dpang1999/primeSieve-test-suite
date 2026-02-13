package finitefft

import (
	"algos/helpers"
	"fmt"
)

// FFTFloat64 is a specialized FFT for real/imag arrays (float64)
type FFT struct{}

// Modular inverse using Extended Euclidean Algorithm
func modInverse(a, p int) int {
	var t, newT int = 0, 1
	var r, newR int = p, a % p

	for newR != 0 {
		quotient := r / newR
		t, newT = newT, t-quotient*newT
		r, newR = newR, r-quotient*newR
	}

	if r > 1 {
		panic(fmt.Sprintf("No modular inverse exists for %d mod %d", a, p))
	}
	if t < 0 {
		t += p
	}
	return t
}

func modPow(base, exp, modulus int) int {
	if modulus == 0 {
		panic("Modulus must be positive")
	}

	result := 1
	base %= modulus

	for exp > 0 {
		if exp%2 == 1 {
			result = (result * base) % modulus
		}
		base = (base * base) % modulus
		exp /= 2
	}

	return result
}

func primitiveRoots(modulus int) int {
	factors := factorize(modulus - 1)
	for g := 2; g < modulus; g++ {
		isRoot := true
		for _, factor := range factors {
			if modPow(g, (modulus-1)/factor, modulus) == 1 {
				isRoot = false
				break
			}
		}
		if isRoot {
			return g
		}
	}

	return 0 // No primitive root found
}

func precomputeRootsOfUnity(n int, direction int, modulus int) []int {
	if (modulus-1)%n != 0 {
		panic("n must divide p-1 for roots of unity to exist in IntModP")
	}

	g := primitiveRoots(modulus)
	//fmt.Printf("Primitive root: %d\n", g)

	omega := modPow(g, (modulus-1)/n, modulus)
	//fmt.Printf("omega %d\n", omega)
	roots := make([]int, n)
	for k := 0; k < n; k++ {
		exponent := (k*direction + modulus - 1) % (modulus - 1)
		roots[k] = int(modPow(omega, exponent, modulus))
		//fmt.Printf("Root %d: %d, exponent: %d\n", k, roots[k], exponent)
	}
	return roots
}

// Factorize a number into its prime factors
func factorize(n int) []int {
	factors := []int{}
	for i := 2; i*i <= n; i++ {
		for n%i == 0 {
			factors = append(factors, i)
			n /= i
		}
	}
	if n > 1 {
		factors = append(factors, n)
	}
	return factors
}

// Transform performs the FFT on the given data (in-place)
// real and imag are both length n (n must be a power of 2)
func (fft FFT) Transform(data []int64, modulus int) {
	fft.transformInternal(data, -1, modulus)
}

// Inverse performs the inverse FFT on the given data (in-place)
func (fft FFT) Inverse(data []int64, modulus int) {
	fft.transformInternal(data, 1, modulus)

	nd := len(data)
	n := nd

	for i := 0; i < nd; i++ {
		data[i] = data[i] * int64(modInverse(n, modulus)) % int64(modulus)
	}
}

// Test performs a round-trip FFT and inverse FFT, returning the RMS error
func (fft FFT) Test(data []int64, modulus int) int {
	n := len(data)
	copyData := make([]int64, n)
	copy(copyData, data)

	fft.Transform(data, modulus)
	fmt.Printf("After transform: %v\n", data)

	fft.Inverse(data, modulus)
	fmt.Printf("After inverse: %v\n", data)

	return 0
}

// transformInternal performs the FFT or inverse FFT
func (fft FFT) transformInternal(data []int64, direction int, modulus int) {
	n := len(data)
	if n <= 1 {
		return
	}

	logn := log2(n)

	bitreverse(data)
	// print bitreverse
	//fmt.Printf("After bit-reverse: %v\n", data)

	roots := precomputeRootsOfUnity(n, direction, modulus)
	// print roots
	//fmt.Printf("Roots of unity: %v\n", roots)

	dual := 1
	for bit := 0; bit < logn; bit++ {
		for a := 0; a < dual; a++ {
			w := int64(roots[a*(n/(2*dual))])
			for b := 0; b < n; b += 2 * dual {
				i := b + a
				j := b + a + dual

				wd := w * data[j] % int64(modulus)
				data[j] = (data[i] + int64(modulus) - wd) % int64(modulus)
				data[i] = (data[i] + wd) % int64(modulus)
			}
		}
		dual *= 2
	}
}

// log2 computes the base-2 logarithm of n, ensuring n is a power of 2
func log2(n int) int {
	log := 0
	k := 1
	for k < n {
		k *= 2
		log++
	}
	if n != (1 << log) {
		panic(fmt.Sprintf("FFT: Data length is not a power of 2!: %d", n))
	}
	return log
}

// bitreverseFloat64 reorders the real/imag arrays in bit-reversed order
func bitreverse(data []int64) {
	n := len(data)
	nm1 := n - 1
	i, j := 0, 0
	for i < nm1 {
		if i < j {
			data[i], data[j] = data[j], data[i]

		}
		k := n >> 1
		for k <= j {
			j -= k
			k >>= 1
		}
		j += k
		i++
	}
}

func TestFFT(n int) {
	// let mode = 0 be for testing
	mode := 1

	if n <= 0 {
		n = 16 // default size
	}
	if mode == 0 {
		rand := helpers.NewLCG(12345, 1345, 16645, 1013904)
		//rand := helpers.NewLCG(12345, 1345, 16645, 1013904)
		data1 := make([]int64, n)
		for i := 0; i < n; i++ {
			data1[i] = int64(rand.NextInt())
		}
		fft := FFT{}
		modulus := helpers.FindCongruentPrime(n)

		fmt.Printf("Go Specialized f64 FFT Test: n=%d\n", n)
		for i := 0; i < 10; i++ {
			fft.transformInternal(data1, -1, modulus)
			fft.transformInternal(data1, 1, modulus)
			fmt.Printf("loop %d done\n", i)
		}
	} else {
		in1 := []int64{38, 0, 44, 87, 6, 45, 22, 93, 0, 0, 0, 0, 0, 0, 0, 0}
		in2 := []int64{80, 18, 62, 90, 17, 96, 27, 97, 0, 0, 0, 0, 0, 0, 0, 0}
		//out := []int64{3040, 684, 5876, 11172, 5420, 16710, 12546, 20555, 16730, 15704, 21665, 5490, 13887, 4645, 9021, 0}
		prime := 40961

		fft := FFT{}
		fft.Transform(in1, prime)
		fft.Transform(in2, prime)
		fmt.Println("Transformed in1: ", in1)
		fmt.Println("Transformed in2: ", in2)

		product := make([]int64, len(in1))

		for i := 0; i < len(in1); i++ {
			product[i] = in1[i] * in2[i] % int64(prime)
		}
		fmt.Println("Product: ", product)

		fft.Inverse(product, prime)
		fmt.Println("Inverse product: ", product)
	}
}
