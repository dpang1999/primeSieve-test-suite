package generic

type IMath[T any] interface {
	abs() T
	sqrt() T
}
