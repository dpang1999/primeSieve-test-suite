package generic

type IPrimitiveRoots[T any] interface {
	PrimitiveRoots(n int64) T
	pow(exp int64) T
}
