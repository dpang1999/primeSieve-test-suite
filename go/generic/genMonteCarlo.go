package generic

import (
	"algos/helpers"
	"fmt"
	"math"
)

// integrate performs the Monte Carlo integration
func integrate[T interface {
	IField[T]
	IOrdered[T]
}](t T, numSamples int) float64 {
	rand := helpers.NewLCG(12345, 1345, 16645, 1013904)
	underCurve := 0
	for count := 0; count < numSamples; count++ {
		x := t.coerceFromFloat(rand.NextDouble())
		y := t.coerceFromFloat(rand.NextDouble())

		// Check if x^2 + y^2 <= 1
		if x.m(x).a(y.m(y)).le(t.one()) {
			underCurve++
		}
	}

	// Return the result: (underCurve / numSamples) * 4
	return float64(underCurve) / float64(numSamples) * 4.0
}

func TestGenMonteCarlo(mode int, numSamples int) {

	if mode == 0 {
		// Use DoubleField
		t := DoubleField{Value: 0.0}
		pi := integrate(t, numSamples)
		fmt.Printf("Go generic doublefield monte carlo\n")
		fmt.Printf("Pi is approximately: %.6f\n", pi)
		fmt.Printf("Num samples: %d\n", numSamples)
		fmt.Printf("RMS Error: %.6f\n", math.Abs(math.Pi-pi))
	} else if mode == 1 {
		// Use SingleField
		t := SingleField{Value: 0}
		pi := integrate(t, numSamples)
		fmt.Printf("Go generic singlefield monte carlo\n")
		fmt.Printf("Pi is approximately: %.6f\n", pi)
		fmt.Printf("Num samples: %d\n", numSamples)
		fmt.Printf("RMS Error: %.6f\n", math.Abs(math.Pi-pi))
	} else {
		t := IntModP{Value: 0}
		SetModulus(1000003) // Set a large prime modulus
		pi := integrate(t, numSamples)
		fmt.Printf("Go generic intmodp monte carlo\n")
		fmt.Printf("Pi is approximately: %.6f\n", pi)
		fmt.Printf("Num samples: %d\n", numSamples)
		fmt.Printf("RMS Error: %.6f\n", math.Abs(math.Pi-pi))
	}
}
