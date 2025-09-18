package generic

import (
    "fmt"
    "math"
    "math/cmplx"
    "math/rand"
)

// GenFFT represents the FFT implementation
type GenFFT[N IField[N] & IMath[N] & IPrimitiveRoots[N]] struct {
    c N
}

// NewGenFFT creates a new GenFFT instance
func NewGenFFT[N IField[N] & IMath[N] & IPrimitiveRoots[N]](data N) *GenFFT[N] {
    return &GenFFT[N]{c: data}
}

// Transform performs the FFT on the given data
func (fft *GenFFT[N]) Transform(data []N) {
    fft.transformInternal(data, -1)
}

// Inverse performs the inverse FFT on the given data
func (fft *GenFFT[N]) Inverse(data []N) {
    fft.transformInternal(data, 1)
    n := len(data)
    norm := fft.c.coerceFromFloat(float64(n))
    for i := range data {
        data[i] = data[i].d(norm)
    }
}

// Test performs a round-trip FFT and inverse FFT, returning the RMS error
func (fft *GenFFT[N]) Test(data []N) float64 {
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
        data[i] = NewComplexField(
            fft.c.coerceFromFloat(rand.Float64()),
            fft.c.coerceFromFloat(rand.Float64()),
        )
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

    roots := fft.c.precomputeRootsOfUnity(uint32(n), int32(direction))

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