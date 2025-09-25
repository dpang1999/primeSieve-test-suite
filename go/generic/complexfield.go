package generic

import (
	"fmt"
	"math"
)

// ComplexField represents a complex number with real and imaginary parts of type T
type ComplexField[T interface {
	IField[T]
	IOrdered[T]
	IMath[T]
}] struct {
	Re T
	Im T
}

func (c ComplexField[T]) String() string {
	return fmt.Sprintf("(%.6f + %.6fi)", c.Re.coerceToFloat(), c.Im.coerceToFloat())
}

func (c ComplexField[T]) a(o ComplexField[T]) ComplexField[T] {
	return ComplexField[T]{Re: c.Re.a(o.Re), Im: c.Im.a(o.Im)}
}

func (c ComplexField[T]) s(o ComplexField[T]) ComplexField[T] {
	return ComplexField[T]{Re: c.Re.s(o.Re), Im: c.Im.s(o.Im)}
}

func (c ComplexField[T]) m(o ComplexField[T]) ComplexField[T] {
	// (a + bi)(c + di) = (ac - bd) + (ad + bc)i
	ac := c.Re.m(o.Re)
	bd := c.Im.m(o.Im)
	ad := c.Re.m(o.Im)
	bc := c.Im.m(o.Re)
	return ComplexField[T]{Re: ac.s(bd), Im: ad.a(bc)}
}

func (c ComplexField[T]) d(o ComplexField[T]) ComplexField[T] {
	// (a + bi) / (c + di) = [(ac + bd) + (bc - ad)i] / (c^2 + d^2)
	denom := o.Re.m(o.Re).a(o.Im.m(o.Im))
	if denom.isZero() {
		panic("division by zero")
	}
	ac := c.Re.m(o.Re)
	bd := c.Im.m(o.Im)
	bc := c.Im.m(o.Re)
	ad := c.Re.m(o.Im)
	return ComplexField[T]{
		Re: ac.a(bd).d(denom),
		Im: bc.s(ad).d(denom),
	}
}

func (c ComplexField[T]) coerceFromInt(i int) ComplexField[T] {
	return ComplexField[T]{Re: c.Re.coerceFromInt(i), Im: c.Im.zero()}
}

func (c ComplexField[T]) coerceFromFloat(f float64) ComplexField[T] {
	return ComplexField[T]{Re: c.Re.coerceFromFloat(f), Im: c.Im.zero()}
}

func (c ComplexField[T]) coerceToFloat() float64 {
	return c.Re.coerceToFloat() // Return the real part as a float
}

func (c ComplexField[T]) isZero() bool {
	return c.Re.isZero() && c.Im.isZero()
}

func (c ComplexField[T]) isOne() bool {
	return c.Re.isOne() && c.Im.isZero()
}

func (c ComplexField[T]) zero() ComplexField[T] {
	return ComplexField[T]{Re: c.Re.zero(), Im: c.Im.zero()}
}

func (c ComplexField[T]) one() ComplexField[T] {
	return ComplexField[T]{Re: c.Re.one(), Im: c.Im.zero()}
}

func (c ComplexField[T]) abs() ComplexField[T] {
	magnitude := c.Re.m(c.Re).a(c.Im.m(c.Im)).sqrt().abs()
	return ComplexField[T]{Re: magnitude, Im: c.Im.zero()}
}

func (c ComplexField[T]) sqrt() ComplexField[T] {
	// panic
	panic("Square root not implemented for ComplexField")
}

func (c ComplexField[T]) primitiveRoots(n int64) ComplexField[T] {
	if n == 0 {
		panic("n must be positive")
	}
	angle := 2.0 * math.Pi / float64(n)
	real := c.Re.coerceFromFloat(math.Cos(angle))
	imag := c.Im.coerceFromFloat(math.Sin(angle))
	return ComplexField[T]{Re: real, Im: imag}
}

func (c ComplexField[T]) pow(exponent int64) ComplexField[T] {
	if exponent == 0 {
		return c.one()
	}

	r := c.abs().Re.coerceToFloat()                                 // Modulus
	theta := math.Atan2(c.Im.coerceToFloat(), c.Re.coerceToFloat()) // Argument (angle)

	// Compute new modulus and argument
	newR := math.Pow(r, float64(exponent))
	newTheta := theta * float64(exponent)

	// Convert back to rectangular form
	real := c.Re.coerceFromFloat(newR * math.Cos(newTheta))
	imag := c.Im.coerceFromFloat(newR * math.Sin(newTheta))
	return ComplexField[T]{Re: real, Im: imag}
}

// PrecomputeRootsOfUnity precomputes the roots of unity
func (c ComplexField[T]) precomputeRootsOfUnity(n int, direction int) []ComplexField[T] {
	roots := make([]ComplexField[T], n)
	for k := 0; k < n; k++ {
		angle := 2.0 * math.Pi * float64(k) / float64(n) * float64(direction)
		real := c.Re.coerceFromFloat(math.Cos(angle))
		imag := c.Im.coerceFromFloat(math.Sin(angle))
		roots[k] = ComplexField[T]{Re: real, Im: imag}
	}
	return roots
}

// Implement IOrdered interface for ComplexField
func (c ComplexField[T]) lt(o ComplexField[T]) bool {
	if c.Re.lt(o.Re) {
		return true
	}
	if c.Re.eq(o.Re) {
		return c.Im.lt(o.Im)
	}
	return false
}

func (c ComplexField[T]) le(o ComplexField[T]) bool {
	return c.lt(o) || c.eq(o)
}

func (c ComplexField[T]) gt(o ComplexField[T]) bool {
	if c.Re.gt(o.Re) {
		return true
	}
	if c.Re.eq(o.Re) {
		return c.Im.gt(o.Im)
	}
	return false
}

func (c ComplexField[T]) ge(o ComplexField[T]) bool {
	return c.gt(o) || c.eq(o)
}

func (c ComplexField[T]) eq(o ComplexField[T]) bool {
	return c.Re.eq(o.Re) && c.Im.eq(o.Im)
}

func (c ComplexField[T]) copy() ComplexField[T] {
	return ComplexField[T]{
		Re: c.Re,
		Im: c.Im,
	}
}
