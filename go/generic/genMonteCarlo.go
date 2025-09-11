package generic

import (
	"fmt"
	"math"
	"math/rand"
)

// integrate performs the Monte Carlo integration
func integrate[T IField[T]](t T, numSamples int) T {

	underCurve := 0
	for count := 0; count < numSamples; count++ {
		x := t.coerceFromFloat(rand.Float64())
		y := t.coerceFromFloat(rand.Float64())

		// Check if x^2 + y^2 <= 1
		if x.m(x).a(y.m(y)).coerceToFloat() <= 1.0 {
			underCurve++
		}
	}

	// Return the result: (underCurve / numSamples) * 4
	return t.coerceFromInt(underCurve).
		d(t.coerceFromInt(numSamples)).
		m(t.coerceFromFloat(4.0))
}

func Main2() {
	mode := 0
	numSamples := 1000000

	// Parse command-line arguments (if any)
	// For simplicity, this example assumes no arguments are passed.

	if mode == 0 {
		// Use DoubleField
		t := DoubleField{Value: 0.0}
		pi := integrate(t, numSamples)
		fmt.Printf("Pi is approximately: %.6f\n", pi.coerceToFloat())
		fmt.Printf("Num samples: %d\n", numSamples)
		fmt.Printf("RMS Error: %.6f\n", math.Abs(math.Pi-pi.coerceToFloat()))
	} else {
		// Placeholder for other modes (e.g., IntModP)
		fmt.Println("Other modes are not implemented in this example.")
	}
}
