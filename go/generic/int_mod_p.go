package generic

import (
	"fmt"
)

// IntModP represents an integer modulo a prime number
type IntModP struct {
	Value      uint64
	Modulus    uint64
	PrintShort bool
}

// NewIntModP creates a new IntModP instance
func NewIntModP(value, modulus uint64) IntModP {
	return IntModP{
		Value:      value % modulus,
		Modulus:    modulus,
		PrintShort: true,
	}
}

// Modular inverse using Extended Euclidean Algorithm
func modInverse(a, p uint64) uint64 {
	var t, newT int64 = 0, 1
	var r, newR int64 = int64(p), int64(a % p)

	for newR != 0 {
		quotient := r / newR
		t, newT = newT, t-quotient*newT
		r, newR = newR, r-quotient*newR
	}

	if r > 1 {
		panic(fmt.Sprintf("No modular inverse exists for %d mod %d", a, p))
	}
	if t < 0 {
		t += int64(p)
	}
	return uint64(t)
}

// Implement IField interface for IntModP
func (i IntModP) a(o IntModP) IntModP {
	return NewIntModP(i.Value+o.Value, i.Modulus)
}

func (i *IntModP) ae(o IntModP) {
	i.Value = (i.Value + o.Value) % i.Modulus
}

func (i IntModP) s(o IntModP) IntModP {
	return NewIntModP(i.Value+i.Modulus-o.Value, i.Modulus)
}

func (i *IntModP) se(o IntModP) {
	i.Value = (i.Value + i.Modulus - o.Value) % i.Modulus
}

func (i IntModP) m(o IntModP) IntModP {
	return NewIntModP(i.Value*o.Value, i.Modulus)
}

func (i *IntModP) me(o IntModP) {
	i.Value = (i.Value * o.Value) % i.Modulus
}

func (i IntModP) d(o IntModP) IntModP {
	if o.Value == 0 {
		panic("Division by zero in IntModP")
	}
	inv := modInverse(o.Value, i.Modulus)
	return NewIntModP(i.Value*inv, i.Modulus)
}

func (i *IntModP) de(o IntModP) {
	if o.Value == 0 {
		panic("Division by zero in IntModP")
	}
	inv := modInverse(o.Value, i.Modulus)
	i.Value = (i.Value * inv) % i.Modulus
}

func (i IntModP) coerceFromInt(v int) IntModP {
	return NewIntModP(uint64(v), i.Modulus)
}

func (i IntModP) coerceFromFloat(f float64) IntModP {
	return NewIntModP(uint64(f), i.Modulus)
}

func (i IntModP) coerceToFloat() float64 {
	return float64(i.Value)
}

func (i IntModP) isZero() bool {
	return i.Value == 0
}

func (i IntModP) isOne() bool {
	return i.Value == 1
}

func (i IntModP) zero() IntModP {
	return NewIntModP(0, i.Modulus)
}

func (i IntModP) one() IntModP {
	return NewIntModP(1, i.Modulus)
}

// Implement IMath interface for IntModP
func (i IntModP) abs() float64 {
	return float64(i.Value)
}

func (i *IntModP) sqrt() {
	panic("Square root not implemented for IntModP")
}

// Factorize a number into its prime factors
func factorize(n uint64) []uint64 {
	factors := []uint64{}
	for i := uint64(2); i*i <= n; i++ {
		for n%i == 0 {
			factors = append(factors, i)
			n /= i
		}
	}
	if n > 1 {
		factors = append(factors, n)
	}
	return factors
}

// Modular exponentiation
func modPow(base, exp, modulus uint64) uint64 {
	if modulus == 0 {
		panic("Modulus must be positive")
	}

	result := uint64(1)
	base %= modulus

	for exp > 0 {
		if exp%2 == 1 {
			result = (result * base) % modulus
		}
		base = (base * base) % modulus
		exp /= 2
	}

	return result
}

// Implement IPrimitiveRoots interface for IntModP
func (i IntModP) primitiveRoot(n uint64) IntModP {
	if n == 0 || n >= i.Modulus {
		panic("n must be in range [1, p-1]")
	}

	factors := factorize(i.Modulus - 1)

	for g := uint64(2); g < i.Modulus; g++ {
		isRoot := true
		for _, factor := range factors {
			if modPow(g, (i.Modulus-1)/factor, i.Modulus) == 1 {
				isRoot = false
				break
			}
		}
		if isRoot {
			return NewIntModP(g, i.Modulus)
		}
	}

	return NewIntModP(0, i.Modulus) // No primitive root found
}

func (i IntModP) pow(exp uint64) IntModP {
	return NewIntModP(modPow(i.Value, exp, i.Modulus), i.Modulus)
}

func (i IntModP) precomputeRootsOfUnity(n uint32, direction int32) []IntModP {
	if (i.Modulus-1)%uint64(n) != 0 {
		panic("n must divide p-1 for roots of unity to exist in IntModP")
	}

	g := i.primitiveRoot(i.Modulus - 1)
	omega := g.pow((i.Modulus - 1) / uint64(n))

	roots := make([]IntModP, n)
	for k := uint32(0); k < n; k++ {
		exponent := uint64(k) * uint64(direction) % (i.Modulus - 1)
		roots[k] = omega.pow(exponent)
	}

	return roots
}

// Implement fmt.Stringer interface for IntModP
func (i IntModP) String() string {
	if i.PrintShort {
		return fmt.Sprintf("%d (mod %d)", i.Value, i.Modulus)
	}
	return fmt.Sprintf("IntModP(%d, %d)", i.Value, i.Modulus)
}
