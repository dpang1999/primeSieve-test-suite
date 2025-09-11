package generic

type IOrdered[T any] interface {
	lt(o T) bool
	le(o T) bool
	gt(o T) bool
	ge(o T) bool
	eq(o T) bool
}
