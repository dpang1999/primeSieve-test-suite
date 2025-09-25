package generic

type IPrimitiveRoots[T any] interface {
	primitiveRoots(n int64) T
	pow(exp int64) T
	precomputeRootsOfUnity(n int, direction int) []T
}
