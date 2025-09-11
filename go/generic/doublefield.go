package generic

import (
	"fmt"
	"math"
)

type DoubleField struct {
	Value float64
}

// Implement the Stringer interface for DoubleField
func (d DoubleField) String() string {
	return fmt.Sprintf("%.6f", d.Value) // Format as a float with 6 decimal places
}

// Implement the IField interface for DoubleField

func (d DoubleField) a(o DoubleField) DoubleField {
	return DoubleField{Value: d.Value + o.Value}
}

func (d DoubleField) ae(o DoubleField) {
	d.Value += o.Value
}

func (d DoubleField) s(o DoubleField) DoubleField {
	return DoubleField{Value: d.Value - o.Value}
}

func (d DoubleField) se(o DoubleField) {
	d.Value -= o.Value
}

func (d DoubleField) m(o DoubleField) DoubleField {
	return DoubleField{Value: d.Value * o.Value}
}

func (d DoubleField) me(o DoubleField) {
	d.Value *= o.Value
}

func (d DoubleField) d(o DoubleField) DoubleField {
	if o.Value == 0 {
		panic("division by zero")
	}
	return DoubleField{Value: d.Value / o.Value}
}

func (d DoubleField) de(o DoubleField) {
	if o.Value == 0 {
		panic("division by zero")
	}
	d.Value /= o.Value
}

func (d DoubleField) coerceFromInt(i int) DoubleField {
	return DoubleField{Value: float64(i)}
}

func (d DoubleField) coerceFromFloat(f float64) DoubleField {
	return DoubleField{Value: f}
}

func (d DoubleField) coerceToFloat() float64 {
	return d.Value
}

func (d DoubleField) isZero() bool {
	return d.Value == 0
}

func (d DoubleField) isOne() bool {
	return d.Value == 1
}

func (d DoubleField) zero() DoubleField {
	return DoubleField{Value: 0}
}

func (d DoubleField) one() DoubleField {
	return DoubleField{Value: 1}
}

// Implement IMath interface for DoubleField
func (d DoubleField) abs() {
	if d.Value < 0 {
		d.Value = -d.Value
	}
}

func (d DoubleField) sqrt() {
	if d.Value < 0 {
		panic("square root of negative number")
	}
	d.Value = math.Sqrt(d.Value)
}
