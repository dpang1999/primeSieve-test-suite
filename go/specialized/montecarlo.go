package specialized

import (
	"math/rand"
)

// Monte Carlo integration specialized for float64 (DoubleField)
func TestMonteCarlo(numSamples int) float64 {
	underCurve := 0
	for count := 0; count < numSamples; count++ {
		x := rand.Float64()
		y := rand.Float64()
		if x*x+y*y <= 1.0 {
			underCurve++
		}
	}
	return float64(underCurve) / float64(numSamples) * 4.0
}
