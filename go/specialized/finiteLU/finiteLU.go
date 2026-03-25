package finitelu

import (
	"algos/helpers"
	"fmt"
	"math"
)

var modulus = int(math.Pow(2, 19) - 1)

func modInverse(a, m int) int {
	var m0, x0, x1 = m, 0, 1
	var aa, mm = a, m
	if m == 1 {
		return 0
	}
	for aa > 1 {
		q := aa / mm
		t := mm
		mm = aa % mm
		aa = t
		tmp := x0
		x0 = x1 - q*x0
		x1 = tmp
	}
	if x1 < 0 {
		x1 += m0
	}
	return (x1 % m0)
}

// LU performs LU decomposition for a possibly non-square matrix, matching the Rust implementation
func factor(a [][]int, pivot []int) int {
	n := len(a)
	m := len(a[0])
	minMN := int(math.Min(float64(m), float64(n)))

	for j := 0; j < minMN; j++ {
		// Find pivot in column j and test for singularity
		jp := j
		t := int(math.Abs(float64(a[j][j])))
		for i := j + 1; i < m; i++ {
			ab := int(math.Abs(float64(a[i][j])))
			if ab > t {
				jp = i
				t = ab
			}
		}
		pivot[j] = jp

		// If zero pivot, factorization fails
		if a[jp][j] == 0 {
			fmt.Println("Matrix is singular")
			return 1
		}

		// Swap rows j and jp if needed
		if jp != j {
			a[j], a[jp] = a[jp], a[j]
		}

		// Compute elements j+1:M of jth column
		if j < m-1 {
			recp := modInverse(a[j][j], modulus)
			for k := j + 1; k < m; k++ {
				a[k][j] = (a[k][j]*recp + modulus) % modulus
			}
		}

		// Rank-1 update to trailing submatrix
		if j < minMN-1 {
			for ii := j + 1; ii < m; ii++ {
				aii_j := a[ii][j]
				for jj := j + 1; jj < n; jj++ {
					a[ii][jj] = (a[ii][jj] - (aii_j*a[j][jj])%modulus + modulus) % modulus
				}
			}
		}
	}
	return 0
}

// solve solves the system using LU and pivot, modifies b in place
func solve(lu [][]int, pvt []int, b []int) {
	m := len(lu)
	n := len(lu[0])
	ii := 0
	for i := 0; i < m; i++ {
		ip := pvt[i]
		sum := b[ip]
		b[ip] = b[i]
		if ii == 0 {
			for j := ii; j < i; j++ {
				sum = (sum - (lu[i][j]*b[j])%modulus + modulus) % modulus
			}
		} else if sum == 0 {
			ii = i
		}
		b[i] = sum
	}

	for i := n - 1; i >= 0; i-- {
		sum := b[i]
		for j := i + 1; j < n; j++ {
			sum = (sum - (lu[i][j]*b[j])%modulus + modulus) % modulus
		}
		b[i] = (sum * modInverse(lu[i][i], modulus)) % modulus
	}
}

func printMatrix(a [][]int) {
	for _, row := range a {
		for _, val := range row {
			fmt.Printf("%d ", val)
		}
		fmt.Println()
	}
}

func printVector(b []int) {
	for _, val := range b {
		fmt.Printf("%d ", val)
	}
	fmt.Println()
}

func TestLU(n int) {
	rand := helpers.NewLCG(12345, 1345, 16645, 1013904)
	if n <= 0 {
		n = 4 // default size
	}
	A := make([][]int, n)
	for i := range A {
		row_sum := 0
		A[i] = make([]int, n)
		for j := range A[i] {
			val := rand.NextInt() % modulus
			row_sum += val
			A[i][j] = val
		}
		A[i][i] = (row_sum + rand.NextInt() + 1) % modulus
	}
	b := make([]int, n)
	for i := range b {
		b[i] = rand.NextInt() % modulus
	}
	//printMatrix(A)
	fmt.Println("Go specialized finite field LU")
	fmt.Println("Matrix size: ", n)
	for i := 0; i < 10; i++ {
		pivot := make([]int, n)
		A_clone := A
		b_clone := b
		factor(A_clone, pivot)
		solve(A_clone, pivot, b_clone)
		fmt.Println("Iteration ", i, " complete")
	}
}
