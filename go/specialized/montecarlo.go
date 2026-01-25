package specialized

import (
	"algos/helpers"
)

// Monte Carlo integration specialized for float64 (DoubleField)
func TestMonteCarlo(numSamples int) float64 {
	rand := helpers.NewLCG(12345, 1345, 16645, 1013904)
	underCurve := 0
	for count := 0; count < numSamples; count++ {
		x := rand.NextDouble()
		y := rand.NextDouble()
		if x*x+y*y <= 1.0 {
			underCurve++
		}
	}
	return float64(underCurve) / float64(numSamples) * 4.0
}
