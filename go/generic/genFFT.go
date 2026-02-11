package generic

import (
	"algos/helpers"
	"fmt"
	"math"
	"math/rand"
)

// GenFFT represents the FFT implementation
type GenFFT[N interface {
	IField[N]
	IMath[N]
	IPrimitiveRoots[N]
	ICopiable[N]
	IOrdered[N]
}] struct {
	c N
}

// NewGenFFT creates a new GenFFT instance
func NewGenFFT[N interface {
	IField[N]
	IMath[N]
	IPrimitiveRoots[N]
	ICopiable[N]
	IOrdered[N]
}](data N) GenFFT[N] {
	return GenFFT[N]{c: data}
}

// Transform performs the FFT on the given data
func (fft GenFFT[N]) Transform(data []N) {
	fft.transformInternal(data, -1)
}

// Inverse performs the inverse FFT on the given data
func (fft GenFFT[N]) Inverse(data []N) {
	fft.transformInternal(data, 1)
	n := len(data)
	norm := fft.c.coerceFromFloat(float64(n))
	for i := range data {
		data[i] = data[i].d(norm)
	}
}

// Test performs a round-trip FFT and inverse FFT, returning the RMS error
func (fft GenFFT[N]) Test(data []N) float64 {
	n := len(data)
	copyData := make([]N, n)
	for i, v := range data {
		copyData[i] = v.copy()
	}

	fft.Transform(data)
	fmt.Println("After transform:", data)
	fft.Inverse(data)
	fmt.Println("After inverse:", data)

	var diff float64
	for i := 0; i < n; i++ {
		realDiff := data[i].coerceToFloat() - copyData[i].coerceToFloat()
		diff += realDiff * realDiff
	}
	return math.Sqrt(diff / float64(n))
}

// MakeRandom generates random complex data
func (fft *GenFFT[N]) MakeRandom(n int) []ComplexField[N] {
	data := make([]ComplexField[N], n)
	for i := 0; i < n; i++ {
		data[i] = ComplexField[N]{
			Re: fft.c.coerceFromFloat(rand.Float64()),
			Im: fft.c.coerceFromFloat(rand.Float64()),
		}
	}
	return data
}

// transformInternal performs the FFT or inverse FFT
func (fft *GenFFT[N]) transformInternal(data []N, direction int) {
	n := len(data)
	if n <= 1 {
		return
	}
	logn := log2(n)

	bitreverse(data)

	roots := fft.c.precomputeRootsOfUnity(n, direction)
	//fmt.Println("Roots of unity:", roots)

	dual := 1
	for bit := 0; bit < logn; bit++ {
		for a := 0; a < dual; a++ {
			w := roots[a*(n/(2*dual))]
			for b := 0; b < n; b += 2 * dual {
				i := b + a
				j := b + a + dual

				wd := w.m(data[j])
				data[j] = data[i].s(wd)
				data[i] = data[i].a(wd)
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

// bitreverse reorders the data array in bit-reversed order
func bitreverse[N any](data []N) {
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

func TestGenFFT(n int, fieldType int) {
	var mode = 1 // let mode = 0 be for testing
	// fieldType: 1 = ComplexField[DoubleField], 2 = IntModP

	if mode != 0 {
		if n <= 0 {
			n = 16
		}
		rand := helpers.NewLCG(12345, 1345, 16645, 1013904)
		if fieldType == 1 {
			fmt.Println("Go Generic FFT Complex Field, ", n)
			varRandomNumbers := make([]ComplexField[DoubleField], n)
			for i := 0; i < n; i++ {
				varRandomNumbers[i] = ComplexField[DoubleField]{
					Re: DoubleField{Value: rand.NextDouble()},
					Im: DoubleField{Value: rand.NextDouble()},
				}
			}
			fft := NewGenFFT(ComplexField[DoubleField]{Re: DoubleField{Value: 3.0}, Im: DoubleField{Value: 4.0}})
			for j := 0; j < 10; j++ {
				fft.transformInternal(varRandomNumbers, -1)
				fft.transformInternal(varRandomNumbers, 1)
				fmt.Println("loop ", j, " done")
			}
		} else {
			fmt.Println("Go Generic FFT IntModP, ", n)
			varRandomNumbers := make([]IntModP, n)
			for i := 0; i < n; i++ {
				varRandomNumbers[i] = NewIntModP(uint64(rand.NextInt()))
			}
			fft := NewGenFFT(NewIntModP(3))
			for j := 0; j < 10; j++ {
				fft.transformInternal(varRandomNumbers, -1)
				fft.transformInternal(varRandomNumbers, 1)
				fmt.Println("loop ", j, " done")
			}

		}

	} else {
		rand := helpers.NewLCG(12345, 1345, 16645, 1013904)
		var randomNumbers [10]int
		var randomDoubles [10]float64
		for i := 0; i < 10; i++ {
			randomNumbers[i] = rand.NextInt()
			randomDoubles[i] = rand.NextDouble()
		}
		fmt.Println("Random Integers:", randomNumbers)
		fmt.Println("Random Doubles:", randomDoubles)

		c := ComplexField[DoubleField]{Re: DoubleField{Value: 3.0}, Im: DoubleField{Value: 4.0}}
		fft := NewGenFFT(c)
		n = 268435456

		data1 := make([]ComplexField[DoubleField], 0, n)
		data1 = append(data1, ComplexField[DoubleField]{
			Re: DoubleField{Value: 0.3618031071604718},
			Im: DoubleField{Value: 0.932993485288541},
		})
		data1 = append(data1, ComplexField[DoubleField]{
			Re: DoubleField{Value: 0.8330913489710237},
			Im: DoubleField{Value: 0.32647575623792624},
		})
		data1 = append(data1, ComplexField[DoubleField]{
			Re: DoubleField{Value: 0.2355237906476252},
			Im: DoubleField{Value: 0.34911535662488336},
		})
		data1 = append(data1, ComplexField[DoubleField]{
			Re: DoubleField{Value: 0.4480776326931518},
			Im: DoubleField{Value: 0.6381529437838686},
		})
		fmt.Println("Input data:", data1)
		fft.Transform(data1)
		fmt.Println("After FFT:", data1)
		fft.Inverse(data1)
		fmt.Println("After Inverse FFT:", data1)
		fmt.Printf("RMS Error: %.6f\n", fft.Test(data1))

		// Test with IntModP
		/*
			in1 := []uint64{38, 0, 44, 87, 6, 45, 22, 93, 0, 0, 0, 0, 0, 0, 0, 0}
			in2 := []uint64{80, 18, 62, 90, 17, 96, 27, 97, 0, 0, 0, 0, 0, 0, 0, 0}
			//out := []uint64{3040, 684, 5876, 11172, 5420, 16710, 12546, 20555, 16730, 15704, 21665, 5490, 13887, 4645, 9021, 0}
			prime := 40961
		*/

		in1 := []uint64{11400, 28374, 23152, 9576, 29511, 20787, 13067, 14015, 0, 0, 0, 0, 0, 0, 0, 0}
		in2 := []uint64{30268, 20788, 8033, 15446, 26275, 11619, 2494, 7016, 0, 0, 0, 0, 0, 0, 0, 0}
		//out := []uint64{ 345055200, 1095807432, 1382179648, 1175142886, 2016084656, 2555168834,2179032777, 1990011337, 1860865174, 1389799087, 942120918, 778961552,341270975, 126631482, 98329240, 0}
		prime := 3221225473
		SetModulus(uint64(prime))

		finite := NewIntModP(3)
		fftInt := NewGenFFT(finite)

		dataA := make([]IntModP, 0, len(in1))
		dataB := make([]IntModP, 0, len(in2))
		for i := 0; i < len(in1); i += 1 {
			dataA = append(dataA, NewIntModP(in1[i]))
			dataB = append(dataB, NewIntModP(in2[i]))
		}
		fmt.Println("Input data A:", dataA)
		fmt.Println("Input data B:", dataB)

		fftInt.Transform(dataA)
		fftInt.Transform(dataB)
		fmt.Println("After FFT A:", dataA)
		fmt.Println("After FFT B:", dataB)

		for i := range dataA {
			dataA[i] = dataA[i].m(dataB[i])
		}
		fmt.Println("After pointwise multiplication:", dataA)

		fftInt.Inverse(dataA)
		fmt.Println("After Inverse FFT:", dataA)
	}
}
