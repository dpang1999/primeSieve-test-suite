package generic

type IPrimitiveRoots[T any] interface {
	primitiveRoots(n int64) T
	pow(exp int64) T
	precomputeRootsOfUnity(n uint32, direction int32) []T
}
