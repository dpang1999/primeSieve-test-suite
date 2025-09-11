package generic

import (
	"fmt"
	"math"
)

type SingleField struct {
	Value float32
}

// Implement the Stringer interface for SingleField
func (d SingleField) String() string {
	return fmt.Sprintf("%.6f", d.Value) // Format as a float with 6 decimal places
}

// Implement the IField interface for SingleField

func (d SingleField) a(o SingleField) SingleField {
	return SingleField{Value: d.Value + o.Value}
}

func (d SingleField) ae(o SingleField) {
	d.Value += o.Value
}

func (d SingleField) s(o SingleField) SingleField {
	return SingleField{Value: d.Value - o.Value}
}

func (d SingleField) se(o SingleField) {
	d.Value -= o.Value
}

func (d SingleField) m(o SingleField) SingleField {
	return SingleField{Value: d.Value * o.Value}
}

func (d SingleField) me(o SingleField) {
	d.Value *= o.Value
}

func (d SingleField) d(o SingleField) SingleField {
	if o.Value == 0 {
		panic("division by zero")
	}
	return SingleField{Value: d.Value / o.Value}
}

func (d SingleField) de(o SingleField) {
	if o.Value == 0 {
		panic("division by zero")
	}
	d.Value /= o.Value
}

func (d SingleField) coerceFromInt(i int) SingleField {
	return SingleField{Value: float32(i)}
}

func (d SingleField) coerceFromFloat(f float64) SingleField {
	return SingleField{Value: float32(f)}
}

func (d SingleField) coerceToFloat() float64 {
	return float64(d.Value)
}

func (d SingleField) isZero() bool {
	return d.Value == 0
}

func (d SingleField) isOne() bool {
	return d.Value == 1
}

func (d SingleField) zero() SingleField {
	return SingleField{Value: 0}
}

func (d SingleField) one() SingleField {
	return SingleField{Value: 1}
}

// Implement IMath interface for SingleField

func (d SingleField) abs() {
	d.Value = float32(math.Abs(float64(d.Value)))
}
func (d SingleField) sqrt() {
	d.Value = float32(math.Sqrt(float64(d.Value)))
}
