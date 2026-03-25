package specialized

import (
	"algos/helpers"
	"math"
)

// Monte Carlo integration specialized for float64 (DoubleField)
func TestMonteCarlo(numSamples int) {
	rand := helpers.NewLCG(12345, 1345, 16645, 1013904)
	underCurve := 0
	for count := 0; count < numSamples; count++ {
		x := rand.NextDouble()
		y := rand.NextDouble()
		if x*x+y*y <= 1.0 {
			underCurve++
		}
	}
	pi := float64(underCurve) / float64(numSamples) * 4.0
	println("Go specialized doublefield monte carlo")
	println("Pi is approximately: %.6f\n", pi)
	println("Num samples: %d\n", numSamples)
	println("RMS Error: %.6f\n", math.Abs(math.Pi-pi))
}
