package specialized

import (
	"algos/helpers"
	"fmt"
	"math"
)

// FFTFloat64 is a specialized FFT for real/imag arrays (float64)
type FFT struct{}

// Transform performs the FFT on the given data (in-place)
// real and imag are both length n (n must be a power of 2)
func (fft FFT) Transform(data []float64) {
	fft.transformInternal(data, -1)
}

// Inverse performs the inverse FFT on the given data (in-place)
func (fft FFT) Inverse(data []float64) {
	fft.transformInternal(data, 1)
	nd := len(data)
	n := nd / 2
	norm := 1.0 / float64(n)

	for i := 0; i < nd; i++ {
		data[i] *= norm
	}
}

// Test performs a round-trip FFT and inverse FFT, returning the RMS error
func (fft FFT) Test(data []float64) float64 {
	n := len(data)
	copyData := make([]float64, n)
	copy(copyData, data)

	fft.Transform(data)
	fmt.Printf("After transform: %v\n", data)

	fft.Inverse(data)
	fmt.Printf("After inverse: %v\n", data)

	var diff float64
	for i := 0; i < n; i++ {
		d := data[i] - copyData[i]
		diff += d * d
	}
	return math.Sqrt(diff / float64(n))
}

// transformInternal performs the FFT or inverse FFT
func (fft FFT) transformInternal(data []float64, direction int) {
	n := len(data)
	if n <= 1 {
		return
	}
	n /= 2
	logn := log2(n)

	bitreverse(data)
	// print bitreverse
	//fmt.Printf("After bit-reverse: %v\n", data)

	for bit := 0; bit < logn; bit++ {
		dual := 1 << bit
		w_real := 1.0
		w_imag := 0.0

		theta := 2.0 * float64(direction) * math.Pi / (2.0 * float64(dual))
		s := math.Sin(theta)
		t := math.Sin(theta / 2.0)
		s2 := 2.0 * t * t

		// a = 0
		for b := 0; b < n; b += 2 * dual {
			i := 2 * b
			j := 2 * (b + dual)

			wd_real := data[j]
			wd_imag := data[j+1]
			data[j] = data[i] - wd_real
			data[j+1] = data[i+1] - wd_imag
			data[i] += wd_real
			data[i+1] += wd_imag
		}

		// a = 1 .. (dual-1)
		for a := 1; a < dual; a++ {
			// Trigonometric recurrence for w -> exp(i*theta) * w
			tmp_real := w_real - s*w_imag - s2*w_real
			tmp_imag := w_imag + s*w_real - s2*w_imag
			w_real = tmp_real
			w_imag = tmp_imag

			for b := 0; b < n; b += 2 * dual {
				i := 2 * (b + a)
				j := 2 * (b + a + dual)

				z1_real := data[j]
				z1_imag := data[j+1]
				wd_real := w_real*z1_real - w_imag*z1_imag
				wd_imag := w_real*z1_imag + w_imag*z1_real

				data[j] = data[i] - wd_real
				data[j+1] = data[i+1] - wd_imag
				data[i] += wd_real
				data[i+1] += wd_imag
			}
		}
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
func bitreverse(data []float64) {
	n := len(data) / 2
	nm1 := n - 1
	i, j := 0, 0
	for i < nm1 {
		if i < j {
			ii := 2 * i
			jj := 2 * j
			data[ii], data[jj] = data[jj], data[ii]
			data[ii+1], data[jj+1] = data[jj+1], data[ii+1]
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
	if mode != 0 {
		rand := helpers.NewLCG(12345, 1345, 16645, 1013904)
		data1 := make([]float64, 2*n)
		for i := 0; i < n; i++ {
			data1[2*i] = rand.NextDouble()
			data1[2*i+1] = rand.NextDouble()
		}
		fft := FFT{}
		rmsError := fft.Test(data1)
		fmt.Printf("FFT RMS Error: %e\n", rmsError)
	} else {

		data1 := []float64{
			0.3618031071604718, 0.932993485288541,
			0.8330913489710237, 0.32647575623792624,
			0.2355237906476252, 0.34911535662488336,
			0.4480776326931518, 0.6381529437838686,
		}

		fft := FFT{}
		rmsError := fft.Test(data1)
		fmt.Printf("FFT RMS Error: %e\n", rmsError)
	}
}
