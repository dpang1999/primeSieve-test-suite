package generic

import "fmt"

// VecExponents is a simple slice-based exponent vector
// Implements IExponents interface

type VecExponents struct {
	Exponents []uint32
}

func NewVecExponents(exps []uint32) VecExponents {
	return VecExponents{Exponents: append([]uint32{}, exps...)}
}

func (v VecExponents) add(o VecExponents) VecExponents {
	result := make([]uint32, len(v.Exponents))
	for i := range v.Exponents {
		result[i] = v.Exponents[i] + o.Exponents[i]
	}
	return VecExponents{Exponents: result}
}

func (v VecExponents) sub(o VecExponents) VecExponents {
	result := make([]uint32, len(v.Exponents))
	for i := range v.Exponents {
		result[i] = v.Exponents[i] - o.Exponents[i]
	}
	return VecExponents{Exponents: result}
}

func (v VecExponents) lcm(o VecExponents) VecExponents {
	result := make([]uint32, len(v.Exponents))
	for i := range v.Exponents {
		if v.Exponents[i] > o.Exponents[i] {
			result[i] = v.Exponents[i]
		} else {
			result[i] = o.Exponents[i]
		}
	}
	return VecExponents{Exponents: result}
}

func (v VecExponents) deg() int {
	var sum int = 0
	for _, e := range v.Exponents {
		sum += int(e)
	}
	return sum
}

func (v VecExponents) lexCompare(o VecExponents) int {
	for i := range v.Exponents {
		if v.Exponents[i] < o.Exponents[i] {
			return -1
		} else if v.Exponents[i] > o.Exponents[i] {
			return 1
		}
	}
	return 0
}

func (v VecExponents) canReduce(divisor VecExponents) bool {
	for i := range v.Exponents {
		if v.Exponents[i] < divisor.Exponents[i] {
			return false
		}
	}
	return true
}

func (v VecExponents) equals(o VecExponents) bool {
	if len(v.Exponents) != len(o.Exponents) {
		return false
	}
	for i := range v.Exponents {
		if v.Exponents[i] != o.Exponents[i] {
			return false
		}
	}
	return true
}

func (v VecExponents) String() string {
	degree := v.deg()
	s := fmt.Sprintf("Degree: %04X, Exponents: ", degree)
	for _, exp := range v.Exponents {
		s += fmt.Sprintf("%02X ", exp)
	}
	return s
}
