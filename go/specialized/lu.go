package specialized

import (
	"algos/helpers"
	"fmt"
	"math"
)

// LU performs LU decomposition for a possibly non-square matrix, matching the Rust implementation
func factor(a [][]float64, pivot []int) int {
	n := len(a)
	m := len(a[0])
	minMN := int(math.Min(float64(m), float64(n)))

	for j := 0; j < minMN; j++ {
		// Find pivot in column j and test for singularity
		jp := j
		t := math.Abs(a[j][j])
		for i := j + 1; i < m; i++ {
			ab := math.Abs(a[i][j])
			if ab > t {
				jp = i
				t = ab
			}
		}
		pivot[j] = jp

		// If zero pivot, factorization fails
		if a[jp][j] == 0.0 {
			return 1
		}

		// Swap rows j and jp if needed
		if jp != j {
			a[j], a[jp] = a[jp], a[j]
		}

		// Compute elements j+1:M of jth column
		if j < m-1 {
			recp := 1.0 / a[j][j]
			for k := j + 1; k < m; k++ {
				a[k][j] *= recp
			}
		}

		// Rank-1 update to trailing submatrix
		if j < minMN-1 {
			for ii := j + 1; ii < m; ii++ {
				aii_j := a[ii][j]
				for jj := j + 1; jj < n; jj++ {
					a[ii][jj] -= aii_j * a[j][jj]
				}
			}
		}
	}
	return 0
}

// solve solves the system using LU and pivot, modifies b in place
func solve(lu [][]float64, pvt []int, b []float64) {
	m := len(lu)
	n := len(lu[0])
	ii := 0
	for i := 0; i < m; i++ {
		ip := pvt[i]
		sum := b[ip]
		b[ip] = b[i]
		if ii == 0 {
			for j := ii; j < i; j++ {
				sum -= lu[i][j] * b[j]
			}
		} else if sum == 0.0 {
			ii = i
		}
		b[i] = sum
	}

	for i := n - 1; i >= 0; i-- {
		sum := b[i]
		for j := i + 1; j < n; j++ {
			sum -= lu[i][j] * b[j]
		}
		b[i] = sum / lu[i][i]
	}
}

func printMatrix(a [][]float64) {
	for _, row := range a {
		for _, val := range row {
			fmt.Printf("%.3f ", val)
		}
		fmt.Println()
	}
}

func printVector(b []float64) {
	for _, val := range b {
		fmt.Printf("%.3f ", val)
	}
	fmt.Println()
}

func TestLU() {
	rand := helpers.NewLCG(12345, 1345, 65, 17)
	n := 4
	a := make([][]float64, n)
	for i := range a {
		a[i] = make([]float64, n)
		for j := range a[i] {
			a[i][j] = float64(rand.NextDouble() * 1000) // simple deterministic values, replace with random if needed
		}
	}
	b := make([]float64, n)
	for i := range b {
		b[i] = float64(rand.NextDouble() * 1000) // simple deterministic values, replace with random if needed
	}
	pivot := make([]int, n)

	printMatrix(a)
	factor(a, pivot)
	fmt.Println("b:")
	printVector(b)
	solve(a, pivot, b)
	fmt.Println("Solution (x):")
	printVector(b)
}
