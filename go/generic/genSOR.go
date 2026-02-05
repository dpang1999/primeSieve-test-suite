package generic

import (
	"fmt"
)

// numFlops calculates the number of floating-point operations
func numFlops(m, n, numIterations int) float64 {
	md := float64(m)
	nd := float64(n)
	numIterD := float64(numIterations)
	return (md - 1.0) * (nd - 1.0) * numIterD * 6.0
}

// execute performs the SOR algorithm
func execute[R IField[R]](omega R, g [][]R, numIterations int) {
	m := len(g)
	n := len(g[0])

	// Create constants for the algorithm
	four := omega.one()
	for i := 0; i < 2; i++ {
		four = four.a(four) // The dumbest way to make four
	}
	omegaOverFour := omega.d(four)
	oneMinusOmega := omega.one().s(omega)

	mm1 := m - 1
	nm1 := n - 1

	for iter := 0; iter < numIterations; iter++ {
		for i := 1; i < mm1; i++ {
			for j := 1; j < nm1; j++ {
				up := g[i-1][j]
				down := g[i+1][j]
				left := g[i][j-1]
				right := g[i][j+1]
				center := g[i][j]

				neighborSum := up.a(down).a(left).a(right)
				newVal := omegaOverFour.m(neighborSum).a(oneMinusOmega.m(center))
				g[i][j] = newVal
			}
		}
	}
}

// printMatrix prints a matrix
func printMatrix[R IField[R]](a [][]R) {
	for _, row := range a {
		for _, val := range row {
			fmt.Printf("%v ", val)
		}
		fmt.Println()
	}
	fmt.Println()
}

// primeSieve generates a list of prime numbers up to a given limit
func primeSieve(limit int) []int {
	numbers := make([]bool, limit)
	primes := []int{}
	for i := range numbers {
		numbers[i] = true
	}
	numbers[0], numbers[1] = false, false

	for i := 2; i < limit; i++ {
		if numbers[i] {
			primes = append(primes, i)
			for j := i * i; j < limit; j += i {
				numbers[j] = false
			}
		}
	}

	return primes
}

// main function to run the SOR algorithm
func TestGenSOR() {
	m := 10
	n := 10
	mode := 3 // 1: SingleField, 2: DoubleField, 3: IntModP
	numIterations := 100

	fmt.Println("Running SOR Algorithm")

	if mode == 1 {
		fmt.Println("Using SingleField")
		omega := SingleField{Value: 1.5}
		g := make([][]SingleField, m)
		for i := range g {
			g[i] = make([]SingleField, n)
			for j := range g[i] {
				g[i][j] = omega.zero()
			}
		}

		// Set boundary conditions
		for i := 0; i < m; i++ {
			g[i][0] = omega.zero()   // Left edge
			g[i][n-1] = omega.zero() // Right edge
		}
		for j := 0; j < n; j++ {
			g[0][j] = SingleField{Value: 100.0} // Top edge (hot)
			g[m-1][j] = omega.zero()            // Bottom edge (cold)
		}

		fmt.Println("Initial grid:")
		printMatrix(g)

		execute(omega, g, numIterations)

		fmt.Println("\nSteady-state temperature distribution:")
		printMatrix(g)
	} else if mode == 2 {
		fmt.Println("Using DoubleField")
		omega := DoubleField{Value: 1.5}
		g := make([][]DoubleField, m)
		for i := range g {
			g[i] = make([]DoubleField, n)
			for j := range g[i] {
				g[i][j] = omega.zero()
			}
		}

		// Set boundary conditions
		for i := 0; i < m; i++ {
			g[i][0] = omega.zero()   // Left edge
			g[i][n-1] = omega.zero() // Right edge
		}
		for j := 0; j < n; j++ {
			g[0][j] = DoubleField{Value: 100.0} // Top edge (hot)
			g[m-1][j] = omega.zero()            // Bottom edge (cold)
		}

		fmt.Println("Initial grid:")
		printMatrix(g)

		execute(omega, g, numIterations)

		fmt.Println("\nSteady-state temperature distribution:")
		printMatrix(g)
	} else {
		fmt.Println("Using IntModP")
		primes := primeSieve(46340) // Max sqrt of int32
		prime := primes[len(primes)-1]
		SetModulus(uint64(prime))

		omega := IntModP{Value: 3}.d(IntModP{Value: 2})
		g := make([][]IntModP, m)
		for i := range g {
			g[i] = make([]IntModP, n)
			for j := range g[i] {
				g[i][j] = omega.zero()
			}
		}

		// Set boundary conditions
		for i := 0; i < m; i++ {
			g[i][0] = omega.zero()   // Left edge
			g[i][n-1] = omega.zero() // Right edge
		}
		for j := 0; j < n; j++ {
			g[0][j] = IntModP{Value: 100} // Top edge (hot)
			g[m-1][j] = omega.zero()      // Bottom edge (cold)
		}

		/*fmt.Println("Initial grid:")
		printMatrix(g)

		execute(omega, g, numIterations)

		fmt.Println("\nSteady-state temperature distribution:")
		printMatrix(g)*/
	}
}
