package generic

type IExponents[T any] interface {
	add(o T) T
	sub(o T) T
	lcm(o T) T
	deg() int
	lexCompare(o T) int
	canReduce(o T) bool
	equals(o T) bool
}
