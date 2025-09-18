package generic

import (
	"fmt"
	"math"
	"math/rand"
)

// GenLU represents the LU factorization and pivot vector
type GenLU[R IField[R]] struct {
	LU    [][]R
	Pivot []int
}

// Factor performs LU factorization (in place)
func Factor[R interface {
	IField[R]
	IMath[R]
}](A [][]R, pivot []int) int {
	N := len(A)
	M := len(A[0])
	minMN := int(math.Min(float64(M), float64(N)))

	for j := 0; j < minMN; j++ {
		// Find pivot in column j and test for singularity
		jp := j
		t := A[j][j]
		t = t.abs()

		for i := j + 1; i < M; i++ {
			ab := A[i][j]
			ab = ab.abs()
			if ab.coerceToFloat() > t.coerceToFloat() {
				jp = i
				t = ab
			}
		}
		pivot[j] = jp

		// If zero pivot, factorization fails
		if A[jp][j].coerceToFloat() == 0 {
			fmt.Println("Matrix is singular")
			return 1
		}

		// Swap rows j and jp if needed
		if jp != j {
			A[j], A[jp] = A[jp], A[j]
		}

		// Compute elements j+1:M of jth column
		if j < M-1 {
			recp := A[j][j].one().d(A[j][j])
			for k := j + 1; k < M; k++ {
				A[k][j] = A[k][j].m(recp)
			}
		}

		// Rank-1 update to trailing submatrix
		if j < minMN-1 {
			for ii := j + 1; ii < M; ii++ {
				for jj := j + 1; jj < N; jj++ {
					A[ii][jj] = A[ii][jj].s(A[ii][j].m(A[j][jj]))
				}
			}
		}
	}

	return 0
}

// Solve solves a linear system using a prefactored matrix in LU form
func Solve[R interface {
	IField[R]
	IMath[R]
}](LU [][]R, pivot []int, b []R) []R {
	M := len(LU)
	N := len(LU[0])
	x := make([]R, len(b))
	copy(x, b)

	ii := 0
	for i := 0; i < M; i++ {
		ip := pivot[i]
		sum := x[ip]
		x[ip] = x[i]
		if ii == 0 {
			for j := ii; j < i; j++ {
				sum = sum.s(LU[i][j].m(x[j]))
			}
		} else if sum.coerceToFloat() == 0 {
			ii = i
		}
		x[i] = sum
	}

	for i := N - 1; i >= 0; i-- {
		sum := x[i]
		for j := i + 1; j < N; j++ {
			sum = sum.s(LU[i][j].m(x[j]))
		}
		x[i] = sum.d(LU[i][i])
	}

	return x
}

// MultiplyMatrices multiplies a matrix by a vector
func MultiplyMatrices[R interface {
	IField[R]
	IMath[R]
}](A [][]R, b []R) []R {
	M := len(A)
	N := len(A[0])
	product := make([]R, M)

	for i := 0; i < M; i++ {
		sum := A[0][0].zero()
		for j := 0; j < N; j++ {
			sum = sum.a(A[i][j].m(b[j]))
		}
		product[i] = sum
	}

	return product
}

// PrimeSieve generates a list of prime numbers up to a given limit
func PrimeSieve(limit int) []int {
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

// Run executes the LU factorization and solving process
func Run[R interface {
	IField[R]
	IMath[R]
}](A [][]R, b []R, pivot []int) {
	fmt.Println("Matrix A:")
	PrintMatrix(A)

	Acopy := make([][]R, len(A))
	for i := range A {
		Acopy[i] = make([]R, len(A[i]))
		copy(Acopy[i], A[i])
	}

	Factor(A, pivot)
	fmt.Println("Vector b:")
	PrintVector(b)

	x := Solve(A, pivot, b)
	fmt.Println("Solution x:")
	PrintVector(x)

	product := MultiplyMatrices(Acopy, x)
	fmt.Println("Product Ax:")
	PrintVector(product)
}

// PrintMatrix prints a matrix
func PrintMatrix[R IField[R]](A [][]R) {
	for _, row := range A {
		for _, val := range row {
			fmt.Printf("%v ", val)
		}
		fmt.Println()
	}
	fmt.Println()
}

// PrintVector prints a vector
func PrintVector[R IField[R]](b []R) {
	for _, val := range b {
		fmt.Printf("%v ", val)
	}
	fmt.Println()
}

func TestGenLU() {
	n := 4
	var A [][]DoubleField
	var b []DoubleField
	pivot := make([]int, n)

	// Randomly populate A and b
	A = make([][]DoubleField, n)
	for i := range A {
		A[i] = make([]DoubleField, n)
		for j := range A[i] {
			A[i][j] = DoubleField{Value: rand.Float64() * 10.0}
		}
	}
	b = make([]DoubleField, n)
	for i := range b {
		b[i] = DoubleField{Value: rand.Float64() * 10.0}
	}

	Run(A, b, pivot)
}
