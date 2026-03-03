package generic

type IField[T any] interface {
	a(o T) T
	ae(o T)
	s(o T) T
	se(o T)
	m(o T) T
	me(o T)
	d(o T) T
	de(o T)

	coerceFromInt(i int) T
	coerceFromFloat(f float64) T
	coerceToFloat() float64

	isZero() bool
	isOne() bool
	zero() T
	one() T
}
