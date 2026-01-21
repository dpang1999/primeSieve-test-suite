package specialized

// SOR specialized for float64 (DoubleField)
func SOR(omega float64, g [][]float64, numIterations int) {
	m := len(g)
	n := len(g[0])

	omegaOverFour := omega * 0.25
	oneMinusOmega := 1.0 - omega
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
				neighborSum := up + down + left + right
				g[i][j] = oneMinusOmega*center + omegaOverFour*neighborSum
			}
		}
	}
}

func TestSOR(n int, iterations int) {
	if n <= 0 {
		n = 10 // default size
	}
	if iterations <= 0 {
		iterations = 100 // default iterations
	}
	omega := 1.5
	g := make([][]float64, n)
	for i := range g {
		g[i] = make([]float64, n)
		for j := range g[i] {
			if i == 0 {
				g[i][j] = 100
			}
		}
	}

	//print g
	printMatrix(g)

	SOR(omega, g, iterations)

	printMatrix(g)

}
