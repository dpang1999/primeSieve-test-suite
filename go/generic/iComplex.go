package generic

type IComplex[T any] interface {
	fromPolar(r T, theta T) T
}
