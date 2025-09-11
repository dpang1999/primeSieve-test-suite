package generic

type ICopiable[T any] interface {
	copy() T
}
