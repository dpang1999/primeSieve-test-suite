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
		// For the sake of fairness, modulus will always be 2^32 to bypass modulus bias between languages
	}
}

// Next generates the next random number in the sequence.
func (l *LCG) NextInt() int {
	l.lastNumber = (l.a*l.lastNumber + l.c) // modulus is treated as 2^32, 32 bit overflow will automatically wrap around
	return l.lastNumber & 0x7FFFFFFF        // ignore bit sign bit to ensure non-negative output
}

// NextDouble generates the next random number in the sequence as a float64 in [0, 1).
func (l *LCG) NextDouble() float64 {
	return float64(l.NextInt()) / 4294967296.0 // 2^32
}
