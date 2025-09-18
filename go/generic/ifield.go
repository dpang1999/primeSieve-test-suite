package generic

type IField[T any] interface {
	a(o T) T
	s(o T) T
	m(o T) T
	d(o T) T

	coerceFromInt(i int) T
	coerceFromFloat(f float64) T
	coerceToFloat() float64

	isZero() bool
	isOne() bool
	zero() T
	one() T
}
