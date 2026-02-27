package helpers

type LCG struct {
	lastNumber int
	a          int
	c          int
	m          int
}

// NewLCG creates a new LCG with the given parameters.
func NewLCG(seed, modulo, multiplier, increment int) *LCG {
	return &LCG{
		lastNumber: seed,
		a:          multiplier,
		c:          increment,
		m:          modulo,
	}
}

// Next generates the next random number in the sequence.
func (l *LCG) NextInt() int {
	l.lastNumber = (l.a*l.lastNumber + l.c) % l.m
	return l.lastNumber
}

// NextDouble generates the next random number in the sequence as a float64 in [0, 1).
func (l *LCG) NextDouble() float64 {
	return float64(l.NextInt()) / float64(l.m)
}
