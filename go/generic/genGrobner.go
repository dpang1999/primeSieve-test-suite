package generic

import (
	"fmt"
	"hash/fnv"
	"math"
	"sort"
)

type TermOrder int

const (
	Lex TermOrder = iota
	GrLex
	RevLex
)

var termOrder TermOrder = Lex

type Term[N interface {
	IField[N]
	IMath[N]
	ICopiable[N]
	IOrdered[N]
}, E interface {
	IExponents[E]
}] struct {
	Coefficient N
	Exponents   E
}

type Polynomial[N interface {
	IField[N]
	IMath[N]
	ICopiable[N]
	IOrdered[N]
}, E interface {
	IExponents[E]
}] struct {
	Terms []Term[N, E]
}

func (t Term[N, E]) Compare(other Term[N, E]) int {
	switch termOrder {
	case Lex:
		return t.Exponents.lexCompare(other.Exponents)
	case GrLex:
		tDeg := t.Exponents.deg()
		oDeg := other.Exponents.deg()
		if tDeg != oDeg {
			return tDeg - oDeg
		}
		return t.Exponents.lexCompare(other.Exponents)
	case RevLex:
		tDeg := t.Exponents.deg()
		oDeg := other.Exponents.deg()
		if tDeg != oDeg {
			return tDeg - oDeg
		}
		// For RevLex, reverse lexCompare
		return -t.Exponents.lexCompare(other.Exponents)
	default:
		return 0
	}
}

func NewPolynomial[N interface {
	IField[N]
	IMath[N]
	ICopiable[N]
	IOrdered[N]
}, E interface {
	IExponents[E]
}](terms []Term[N, E]) Polynomial[N, E] {
	filtered := make([]Term[N, E], 0, len(terms))
	for _, t := range terms {
		if math.Abs(t.Coefficient.coerceToFloat()) > 0.0 {
			t.Coefficient = t.Coefficient.coerceFromFloat(
				math.Round(t.Coefficient.coerceToFloat()*1e5) / 1e5)
			filtered = append(filtered, t)
		}
	}
	sort.Slice(filtered, func(i, j int) bool {
		return filtered[i].Compare(filtered[j]) > 0
	})
	return Polynomial[N, E]{Terms: filtered}
}

func (p Polynomial[N, E]) Add(other Polynomial[N, E]) Polynomial[N, E] {
	result := append([]Term[N, E]{}, p.Terms...)
	for _, t := range other.Terms {
		found := false
		for i := range result {
			if result[i].Exponents.equals(t.Exponents) {
				result[i].Coefficient = result[i].Coefficient.a(t.Coefficient)
				found = true
				break
			}
		}
		if !found {
			result = append(result, t)
		}
	}
	return NewPolynomial(result)
}

func (p Polynomial[N, E]) Subtract(other Polynomial[N, E]) Polynomial[N, E] {
	result := append([]Term[N, E]{}, p.Terms...)
	for _, t := range other.Terms {
		found := false
		for i := range result {
			if result[i].Exponents.equals(t.Exponents) {
				result[i].Coefficient = result[i].Coefficient.s(t.Coefficient)
				found = true
				break
			}
		}
		if !found {
			t2 := t
			t2.Coefficient = t2.Coefficient.zero().s(t2.Coefficient)
			result = append(result, t2)
		}
	}
	return NewPolynomial(result)
}

func (p Polynomial[N, E]) MakeMonic() Polynomial[N, E] {
	if len(p.Terms) == 0 {
		return p
	}

	leadCoeff := p.Terms[0].Coefficient
	if leadCoeff.coerceToFloat() == 0 {
		return p
	}

	terms := make([]Term[N, E], len(p.Terms))
	for i, t := range p.Terms {
		terms[i] = Term[N, E]{
			Coefficient: t.Coefficient.d(leadCoeff),
			Exponents:   t.Exponents,
		}
	}

	return NewPolynomial(terms)
}

func (p Polynomial[N, E]) MultiplyByTerm(term Term[N, E]) Polynomial[N, E] {
	terms := make([]Term[N, E], len(p.Terms))
	for i, t := range p.Terms {
		terms[i] = Term[N, E]{
			Coefficient: t.Coefficient.m(term.Coefficient),
			Exponents:   t.Exponents.add(term.Exponents),
		}
	}
	return NewPolynomial(terms)
}

func polynomialHash[N interface {
	IField[N]
	IMath[N]
	ICopiable[N]
	IOrdered[N]
	toBytes() []byte
}, E interface {
	IExponents[E]
	toBytes() []byte
}](terms []Term[N, E]) uint64 {
	h := fnv.New64a()
	for _, term := range terms {
		h.Write(term.Coefficient.toBytes())
		h.Write(term.Exponents.toBytes())
	}
	return h.Sum64()
}

func (p Polynomial[N, E]) Reduce(divisors []Polynomial[N, E]) Polynomial[N, E] {
	result := p
	remainder := []Term[N, E]{}
	for {
		reduced := false
		if len(result.Terms) == 0 {
			break
		}
		lead := result.Terms[0]
		for _, divisor := range divisors {
			if len(divisor.Terms) == 0 {
				continue
			}
			divLead := divisor.Terms[0]
			if lead.Exponents.canReduce(divLead.Exponents) {
				coeff := lead.Coefficient.d(divLead.Coefficient)
				exps := lead.Exponents.sub(divLead.Exponents)
				reductionTerm := Term[N, E]{Coefficient: coeff, Exponents: exps}
				scaledDivisor := divisor.MultiplyByTerm(reductionTerm)
				result = result.Subtract(scaledDivisor)
				reduced = true
				break
			}
		}
		if !reduced {
			remainder = append(remainder, result.Terms[0])
			result.Terms = result.Terms[1:]
		}
	}
	result.Terms = append(result.Terms, remainder...)
	return NewPolynomial(result.Terms)
}

func SPolynomial[N interface {
	IField[N]
	IMath[N]
	ICopiable[N]
	IOrdered[N]
}, E interface {
	IExponents[E]
}](p1 Polynomial[N, E], p2 Polynomial[N, E]) Polynomial[N, E] {
	lead1 := p1.Terms[0]
	lead2 := p2.Terms[0]
	lcmExps := lead1.Exponents.lcm(lead2.Exponents)
	scale1 := lcmExps.sub(lead1.Exponents)
	scale2 := lcmExps.sub(lead2.Exponents)
	scaled1 := p1.MultiplyByTerm(Term[N, E]{Coefficient: lead1.Coefficient.one(), Exponents: scale1})
	scaled2 := p2.MultiplyByTerm(Term[N, E]{Coefficient: lead2.Coefficient.one(), Exponents: scale2})
	result := scaled1.Subtract(scaled2)
	return result
}

func NaiveGrobnerBasis[N interface {
	IField[N]
	IMath[N]
	ICopiable[N]
	IOrdered[N]
}, E interface {
	IExponents[E]
}](polys []Polynomial[N, E]) []Polynomial[N, E] {
	basis := append([]Polynomial[N, E]{}, polys...)
	basisSet := make(map[uint64]struct{})
	for _, poly := range basis {
		key := polynomialHash(poly.Terms)
		basisSet[key] = struct{}{}
	}

	type pair struct {
		i, j int
	}
	var pairs []pair
	for i := 0; i < len(basis); i++ {
		for j := i + 1; j < len(basis); j++ {
			pairs = append(pairs, pair{i, j})
		}
	}
	for len(pairs) > 0 {
		p := pairs[0]
		pairs = pairs[1:]
		spoly := SPolynomial(basis[p.i], basis[p.j])
		reduced := spoly.Reduce(basis)
		if len(reduced.Terms) > 0 {
			key := polynomialHash(reduced.Terms)
			if _, exists := basisSet[key]; !exists {
				basisSet[key] = struct{}{}
				newIdx := len(basis)
				basis = append(basis, reduced)
				for i := 0; i < newIdx; i++ {
					pairs = append(pairs, pair{i, newIdx})
				}
			}
		}
	}

	reducedBasis := []Polynomial[N, E]{}
	for _, poly := range basis {
		basisExcl := make([]Polynomial[N, E], 0, len(basis)-1)
		for _, p := range basis {
			if !polynomialsEqual(p, poly) {
				basisExcl = append(basisExcl, p)
			}
		}
		reduced := poly.Reduce(basisExcl).MakeMonic()
		key := polynomialHash(reduced.Terms)
		if len(reduced.Terms) > 0 {
			found := false
			for _, rb := range reducedBasis {
				if polynomialHash(rb.Terms) == key {
					found = true
					break
				}
			}
			if !found {
				reducedBasis = append(reducedBasis, reduced)
			}
		}
	}
	return reducedBasis
}

func polynomialsEqual[N interface {
	IField[N]
	IMath[N]
	ICopiable[N]
	IOrdered[N]
}, E interface {
	IExponents[E]
}](a, b Polynomial[N, E]) bool {
	if len(a.Terms) != len(b.Terms) {
		return false
	}
	for i := range a.Terms {
		if a.Terms[i].Coefficient.coerceToFloat() != b.Terms[i].Coefficient.coerceToFloat() || a.Terms[i].Exponents.lexCompare(b.Terms[i].Exponents) != 0 {
			return false
		}
	}
	return true
}

/*
func TestGenGrobner[N interface {
	IField[N]
	IMath[N]
	ICopiable[N]
	IOrdered[N]
}, E interface {
	IExponents[E]
}](numPoly int, numTerms int, coeffType int, expType int, orderInt int, modulus int) {
	mode := 1
	if mode != 0 {
		rand := helpers.NewLCG(12345, 1345, 16645, 1013904)
		if numPoly <= 0 {
			numPoly = 3
		}
		if numTerms <= 0 {
			numTerms = 3
		}
		if coeffType < 0 || coeffType > 2 {
			coeffType = 0
		}
		if expType < 0 || expType > 1 {
			expType = 0
		}
		if orderInt < 0 || orderInt > 2 {
			orderInt = 0
		}
		if modulus <= 0 {
			modulus = 13
		}
		order := TermOrder(orderInt)
		numTerms := 3
		var polys []Polynomial[N, E]
		for p := 0; p < numPoly; p++ {
			var terms []Term[N, E]
			for t := 0; t < numTerms; t++ {
				coeff := field.coerceFromFloat(rand.NextDouble()*2 - 1)
				// You must construct exponents of type E here. Example for VecExponents:
				// exps := NewVecExponents([]uint32{...})
				// For BitPackedExponents: exps := NewBitPackedExponents([6]uint8{...})
				// Here, just use expZero as a placeholder (user should replace with actual random exponents)
				terms = append(terms, Term[N, E]{Coefficient: coeff, Exponents: expZero})
			}
			polys = append(polys, NewPolynomial(terms, order))
		}
		basis := NaiveGrobnerBasis(polys, order)
		fmt.Println("Grobner Basis:")
		for i, poly := range basis {
			fmt.Printf("G%d: ", i)
			for _, term := range poly.Terms {
				fmt.Printf("%v*%v ", term.Coefficient, term.Exponents)
			}
			fmt.Println()
		}
	} else {
		fmt.Println("--- Grobner Basis Test ---")
		// Example for VecExponents:
		// p1 := NewPolynomial([]Term[N, E]{
		//     {Coefficient: field.one(), Exponents: NewVecExponents([]uint32{3,0,0})},
		//     ...
		// }, Lex)
		// For BitPackedExponents: use NewBitPackedExponents([6]uint8{...})
		fmt.Println("--- End Grobner Basis Test ---")
	}
}
*/
// GenerateRandomBasis runs Grobner basis computation for all combinations of coefficient and exponent types, similar to the Java main method.
func TestGenGrobner(polyNum, expType, orderInt, mod int) {
	mode := 0 // 0 = test mode with cyclic polynomials, 1 = other testing
	// expType: 0 = vec exponents, 1 = bit-packed exponents
	if mode != 0 {
	} else {
		n := polyNum
		if n <= 0 {
			n = 4
		}
		//modulus = mod
		modulus = uint64(7)
		SetModulus(modulus)
		//termOrder = termOrder
		termOrder = Lex

		if n == 4 {
			if expType == 0 {
				fmt.Println("Go generic finite coeff vec exponent cyclic", n)
				p1 := NewPolynomial([]Term[IntModP, VecExponents]{
					{Coefficient: IntModP{1}, Exponents: NewVecExponents([]uint32{1, 0, 0, 0})},
					{Coefficient: IntModP{1}, Exponents: NewVecExponents([]uint32{0, 1, 0, 0})},
					{Coefficient: IntModP{1}, Exponents: NewVecExponents([]uint32{0, 0, 1, 0})},
					{Coefficient: IntModP{1}, Exponents: NewVecExponents([]uint32{0, 0, 0, 1})},
				})
				p2 := NewPolynomial([]Term[IntModP, VecExponents]{
					{Coefficient: IntModP{1}, Exponents: NewVecExponents([]uint32{1, 1, 0, 0})},
					{Coefficient: IntModP{1}, Exponents: NewVecExponents([]uint32{0, 1, 1, 0})},
					{Coefficient: IntModP{1}, Exponents: NewVecExponents([]uint32{0, 0, 1, 1})},
					{Coefficient: IntModP{1}, Exponents: NewVecExponents([]uint32{1, 0, 0, 1})},
				})
				p3 := NewPolynomial([]Term[IntModP, VecExponents]{
					{Coefficient: IntModP{1}, Exponents: NewVecExponents([]uint32{1, 1, 1, 0})},
					{Coefficient: IntModP{1}, Exponents: NewVecExponents([]uint32{0, 1, 1, 1})},
					{Coefficient: IntModP{1}, Exponents: NewVecExponents([]uint32{1, 0, 1, 1})},
					{Coefficient: IntModP{1}, Exponents: NewVecExponents([]uint32{1, 1, 0, 1})},
				})
				p4 := NewPolynomial([]Term[IntModP, VecExponents]{
					{Coefficient: IntModP{1}, Exponents: NewVecExponents([]uint32{1, 1, 1, 1})},
					{Coefficient: IntModP{modulus - 1}, Exponents: NewVecExponents([]uint32{0, 0, 0, 0})},
				})
				start := []Polynomial[IntModP, VecExponents]{p1, p2, p3, p4}
				for i := 0; i < 10; i++ {
					basis := NaiveGrobnerBasis(start)
					fmt.Printf("Iteration %d: complete\n", i)
					if i == 9 {
						fmt.Println("Final Grobner Basis:")
						for _, poly := range basis {
							fmt.Printf("%v\n\n", poly)
						}
					}
				}
			} else {
				fmt.Println("Go generic finite coeff bitpacked exponent cyclic", n)
				p1 := NewPolynomial([]Term[IntModP, BitPackedExponents]{
					{Coefficient: IntModP{1}, Exponents: NewBitPackedExponents([6]uint8{1, 0, 0, 0, 0, 0})},
					{Coefficient: IntModP{1}, Exponents: NewBitPackedExponents([6]uint8{0, 1, 0, 0, 0, 0})},
					{Coefficient: IntModP{1}, Exponents: NewBitPackedExponents([6]uint8{0, 0, 1, 0, 0, 0})},
					{Coefficient: IntModP{1}, Exponents: NewBitPackedExponents([6]uint8{0, 0, 0, 1, 0, 0})},
				})
				p2 := NewPolynomial([]Term[IntModP, BitPackedExponents]{
					{Coefficient: IntModP{1}, Exponents: NewBitPackedExponents([6]uint8{1, 1, 0, 0, 0, 0})},
					{Coefficient: IntModP{1}, Exponents: NewBitPackedExponents([6]uint8{0, 1, 1, 0, 0, 0})},
					{Coefficient: IntModP{1}, Exponents: NewBitPackedExponents([6]uint8{0, 0, 1, 1, 0, 0})},
					{Coefficient: IntModP{1}, Exponents: NewBitPackedExponents([6]uint8{1, 0, 0, 1, 0, 0})},
				})
				p3 := NewPolynomial([]Term[IntModP, BitPackedExponents]{
					{Coefficient: IntModP{1}, Exponents: NewBitPackedExponents([6]uint8{1, 1, 1, 0, 0, 0})},
					{Coefficient: IntModP{1}, Exponents: NewBitPackedExponents([6]uint8{0, 1, 1, 1, 0, 0})},
					{Coefficient: IntModP{1}, Exponents: NewBitPackedExponents([6]uint8{1, 0, 1, 1, 0, 0})},
					{Coefficient: IntModP{1}, Exponents: NewBitPackedExponents([6]uint8{1, 1, 0, 1, 0, 0})},
				})
				p4 := NewPolynomial([]Term[IntModP, BitPackedExponents]{
					{Coefficient: IntModP{1}, Exponents: NewBitPackedExponents([6]uint8{1, 1, 1, 1, 0, 0})},
					{Coefficient: IntModP{modulus - 1}, Exponents: NewBitPackedExponents([6]uint8{0, 0, 0, 0, 0, 0})},
				})
				start := []Polynomial[IntModP, BitPackedExponents]{p1, p2, p3, p4}
				for i := 0; i < 10; i++ {
					basis := NaiveGrobnerBasis(start)
					fmt.Printf("Iteration %d: complete\n", i)
					if i == 9 {
						fmt.Println("Final Grobner Basis:")
						for _, poly := range basis {
							fmt.Printf("%v\n\n", poly)
						}
					}
				}
			}
		} else if n == 5 {
			if expType == 0 {
				fmt.Println("Go generic finite coeff vec exponent cyclic", n)
				p1 := NewPolynomial([]Term[IntModP, VecExponents]{
					{Coefficient: IntModP{1}, Exponents: NewVecExponents([]uint32{1, 0, 0, 0, 0})},
					{Coefficient: IntModP{1}, Exponents: NewVecExponents([]uint32{0, 1, 0, 0, 0})},
					{Coefficient: IntModP{1}, Exponents: NewVecExponents([]uint32{0, 0, 1, 0, 0})},
					{Coefficient: IntModP{1}, Exponents: NewVecExponents([]uint32{0, 0, 0, 1, 0})},
					{Coefficient: IntModP{1}, Exponents: NewVecExponents([]uint32{0, 0, 0, 0, 1})},
				})
				p2 := NewPolynomial([]Term[IntModP, VecExponents]{
					{Coefficient: IntModP{1}, Exponents: NewVecExponents([]uint32{1, 1, 0, 0, 0})},
					{Coefficient: IntModP{1}, Exponents: NewVecExponents([]uint32{0, 1, 1, 0, 0})},
					{Coefficient: IntModP{1}, Exponents: NewVecExponents([]uint32{0, 0, 1, 1, 0})},
					{Coefficient: IntModP{1}, Exponents: NewVecExponents([]uint32{0, 0, 0, 1, 1})},
					{Coefficient: IntModP{1}, Exponents: NewVecExponents([]uint32{1, 0, 0, 0, 1})},
				})
				p3 := NewPolynomial([]Term[IntModP, VecExponents]{
					{Coefficient: IntModP{1}, Exponents: NewVecExponents([]uint32{1, 1, 1, 0, 0})},
					{Coefficient: IntModP{1}, Exponents: NewVecExponents([]uint32{0, 1, 1, 1, 0})},
					{Coefficient: IntModP{1}, Exponents: NewVecExponents([]uint32{0, 0, 1, 1, 1})},
					{Coefficient: IntModP{1}, Exponents: NewVecExponents([]uint32{1, 0, 0, 1, 1})},
					{Coefficient: IntModP{1}, Exponents: NewVecExponents([]uint32{1, 1, 0, 0, 1})},
				})
				p4 := NewPolynomial([]Term[IntModP, VecExponents]{
					{Coefficient: IntModP{1}, Exponents: NewVecExponents([]uint32{1, 1, 1, 1, 0})},
					{Coefficient: IntModP{1}, Exponents: NewVecExponents([]uint32{0, 1, 1, 1, 1})},
					{Coefficient: IntModP{1}, Exponents: NewVecExponents([]uint32{1, 0, 1, 1, 1})},
					{Coefficient: IntModP{1}, Exponents: NewVecExponents([]uint32{1, 1, 0, 1, 1})},
					{Coefficient: IntModP{1}, Exponents: NewVecExponents([]uint32{1, 1, 1, 0, 1})},
				})
				p5 := NewPolynomial([]Term[IntModP, VecExponents]{
					{Coefficient: IntModP{1}, Exponents: NewVecExponents([]uint32{1, 1, 1, 1, 1})},
					{Coefficient: IntModP{modulus - 1}, Exponents: NewVecExponents([]uint32{0, 0, 0, 0, 0})},
				})
				start := []Polynomial[IntModP, VecExponents]{p1, p2, p3, p4, p5}
				for i := 0; i < 10; i++ {
					basis := NaiveGrobnerBasis(start)
					fmt.Printf("Iteration %d: complete\n", i)
					if i == 9 {
						fmt.Println("Final Grobner Basis:")
						for _, poly := range basis {
							fmt.Printf("%v\n\n", poly)
						}
					}
				}
			} else {
				fmt.Println("Go generic finite coeff bitpacked exponent cyclic", n)
				p1 := NewPolynomial([]Term[IntModP, BitPackedExponents]{
					{Coefficient: IntModP{1}, Exponents: NewBitPackedExponents([6]uint8{1, 0, 0, 0, 0, 0})},
					{Coefficient: IntModP{1}, Exponents: NewBitPackedExponents([6]uint8{0, 1, 0, 0, 0, 0})},
					{Coefficient: IntModP{1}, Exponents: NewBitPackedExponents([6]uint8{0, 0, 1, 0, 0, 0})},
					{Coefficient: IntModP{1}, Exponents: NewBitPackedExponents([6]uint8{0, 0, 0, 1, 0, 0})},
					{Coefficient: IntModP{1}, Exponents: NewBitPackedExponents([6]uint8{0, 0, 0, 0, 1, 0})},
				})
				p2 := NewPolynomial([]Term[IntModP, BitPackedExponents]{
					{Coefficient: IntModP{1}, Exponents: NewBitPackedExponents([6]uint8{1, 1, 0, 0, 0, 0})},
					{Coefficient: IntModP{1}, Exponents: NewBitPackedExponents([6]uint8{0, 1, 1, 0, 0, 0})},
					{Coefficient: IntModP{1}, Exponents: NewBitPackedExponents([6]uint8{0, 0, 1, 1, 0, 0})},
					{Coefficient: IntModP{1}, Exponents: NewBitPackedExponents([6]uint8{0, 0, 0, 1, 1, 0})},
					{Coefficient: IntModP{1}, Exponents: NewBitPackedExponents([6]uint8{1, 0, 0, 0, 1, 0})},
				})
				p3 := NewPolynomial([]Term[IntModP, BitPackedExponents]{
					{Coefficient: IntModP{1}, Exponents: NewBitPackedExponents([6]uint8{1, 1, 1, 0, 0, 0})},
					{Coefficient: IntModP{1}, Exponents: NewBitPackedExponents([6]uint8{0, 1, 1, 1, 0, 0})},
					{Coefficient: IntModP{1}, Exponents: NewBitPackedExponents([6]uint8{0, 0, 1, 1, 1, 0})},
					{Coefficient: IntModP{1}, Exponents: NewBitPackedExponents([6]uint8{1, 0, 0, 1, 1, 0})},
					{Coefficient: IntModP{1}, Exponents: NewBitPackedExponents([6]uint8{1, 1, 0, 0, 1, 0})},
				})
				p4 := NewPolynomial([]Term[IntModP, BitPackedExponents]{
					{Coefficient: IntModP{1}, Exponents: NewBitPackedExponents([6]uint8{1, 1, 1, 1, 0, 0})},
					{Coefficient: IntModP{1}, Exponents: NewBitPackedExponents([6]uint8{0, 1, 1, 1, 1, 0})},
					{Coefficient: IntModP{1}, Exponents: NewBitPackedExponents([6]uint8{1, 0, 1, 1, 1, 0})},
					{Coefficient: IntModP{1}, Exponents: NewBitPackedExponents([6]uint8{1, 1, 0, 1, 1, 0})},
					{Coefficient: IntModP{1}, Exponents: NewBitPackedExponents([6]uint8{1, 1, 1, 0, 1, 0})},
				})
				p5 := NewPolynomial([]Term[IntModP, BitPackedExponents]{
					{Coefficient: IntModP{1}, Exponents: NewBitPackedExponents([6]uint8{1, 1, 1, 1, 1, 0})},
					{Coefficient: IntModP{modulus - 1}, Exponents: NewBitPackedExponents([6]uint8{0, 0, 0, 0, 0, 0})},
				})
				start := []Polynomial[IntModP, BitPackedExponents]{p1, p2, p3, p4, p5}
				for i := 0; i < 10; i++ {
					basis := NaiveGrobnerBasis(start)
					fmt.Printf("Iteration %d: complete\n", i)
					if i == 9 {
						fmt.Println("Final Grobner Basis:")
						for _, poly := range basis {
							fmt.Printf("%v\n\n", poly)
						}
					}
				}
			}
		} else if n == 6 {
			if expType == 0 {
				fmt.Println("Go generic finite coeff vec exponent cyclic", n)
				p1 := NewPolynomial([]Term[IntModP, VecExponents]{
					{Coefficient: IntModP{1}, Exponents: NewVecExponents([]uint32{1, 0, 0, 0, 0, 0})},
					{Coefficient: IntModP{1}, Exponents: NewVecExponents([]uint32{0, 1, 0, 0, 0, 0})},
					{Coefficient: IntModP{1}, Exponents: NewVecExponents([]uint32{0, 0, 1, 0, 0, 0})},
					{Coefficient: IntModP{1}, Exponents: NewVecExponents([]uint32{0, 0, 0, 1, 0, 0})},
					{Coefficient: IntModP{1}, Exponents: NewVecExponents([]uint32{0, 0, 0, 0, 1, 0})},
					{Coefficient: IntModP{1}, Exponents: NewVecExponents([]uint32{0, 0, 0, 0, 0, 1})},
				})
				p2 := NewPolynomial([]Term[IntModP, VecExponents]{
					{Coefficient: IntModP{1}, Exponents: NewVecExponents([]uint32{1, 1, 0, 0, 0, 0})},
					{Coefficient: IntModP{1}, Exponents: NewVecExponents([]uint32{0, 1, 1, 0, 0, 0})},
					{Coefficient: IntModP{1}, Exponents: NewVecExponents([]uint32{0, 0, 1, 1, 0, 0})},
					{Coefficient: IntModP{1}, Exponents: NewVecExponents([]uint32{0, 0, 0, 1, 1, 0})},
					{Coefficient: IntModP{1}, Exponents: NewVecExponents([]uint32{0, 0, 0, 0, 1, 1})},
					{Coefficient: IntModP{1}, Exponents: NewVecExponents([]uint32{1, 0, 0, 0, 0, 1})},
				})
				p3 := NewPolynomial([]Term[IntModP, VecExponents]{
					{Coefficient: IntModP{1}, Exponents: NewVecExponents([]uint32{1, 1, 1, 0, 0, 0})},
					{Coefficient: IntModP{1}, Exponents: NewVecExponents([]uint32{0, 1, 1, 1, 0, 0})},
					{Coefficient: IntModP{1}, Exponents: NewVecExponents([]uint32{0, 0, 1, 1, 1, 0})},
					{Coefficient: IntModP{1}, Exponents: NewVecExponents([]uint32{0, 0, 0, 1, 1, 1})},
					{Coefficient: IntModP{1}, Exponents: NewVecExponents([]uint32{1, 0, 0, 0, 1, 1})},
					{Coefficient: IntModP{1}, Exponents: NewVecExponents([]uint32{1, 1, 0, 0, 0, 1})},
				})
				p4 := NewPolynomial([]Term[IntModP, VecExponents]{
					{Coefficient: IntModP{1}, Exponents: NewVecExponents([]uint32{1, 1, 1, 1, 0, 0})},
					{Coefficient: IntModP{1}, Exponents: NewVecExponents([]uint32{0, 1, 1, 1, 1, 0})},
					{Coefficient: IntModP{1}, Exponents: NewVecExponents([]uint32{0, 0, 1, 1, 1, 1})},
					{Coefficient: IntModP{1}, Exponents: NewVecExponents([]uint32{1, 0, 0, 1, 1, 1})},
					{Coefficient: IntModP{1}, Exponents: NewVecExponents([]uint32{1, 1, 0, 0, 1, 1})},
					{Coefficient: IntModP{1}, Exponents: NewVecExponents([]uint32{1, 1, 1, 0, 0, 1})},
				})
				p5 := NewPolynomial([]Term[IntModP, VecExponents]{
					{Coefficient: IntModP{1}, Exponents: NewVecExponents([]uint32{1, 1, 1, 1, 1, 0})},
					{Coefficient: IntModP{1}, Exponents: NewVecExponents([]uint32{0, 1, 1, 1, 1, 1})},
					{Coefficient: IntModP{1}, Exponents: NewVecExponents([]uint32{1, 0, 1, 1, 1, 1})},
					{Coefficient: IntModP{1}, Exponents: NewVecExponents([]uint32{1, 1, 0, 1, 1, 1})},
					{Coefficient: IntModP{1}, Exponents: NewVecExponents([]uint32{1, 1, 1, 0, 1, 1})},
					{Coefficient: IntModP{1}, Exponents: NewVecExponents([]uint32{1, 1, 1, 1, 0, 1})},
				})
				p6 := NewPolynomial([]Term[IntModP, VecExponents]{
					{Coefficient: IntModP{1}, Exponents: NewVecExponents([]uint32{1, 1, 1, 1, 1, 1})},
					{Coefficient: IntModP{modulus - 1}, Exponents: NewVecExponents([]uint32{0, 0, 0, 0, 0, 0})},
				})
				start := []Polynomial[IntModP, VecExponents]{p1, p2, p3, p4, p5, p6}
				for i := 0; i < 10; i++ {
					basis := NaiveGrobnerBasis(start)
					fmt.Printf("Iteration %d: complete\n", i)
					if i == 9 {
						fmt.Println("Final Grobner Basis:")
						for _, poly := range basis {
							fmt.Printf("%v\n\n", poly)
						}
					}
				}
			} else {
				fmt.Println("Go generic finite coeff bitpacked exponent cyclic", n)
				p1 := NewPolynomial([]Term[IntModP, BitPackedExponents]{
					{Coefficient: IntModP{1}, Exponents: NewBitPackedExponents([6]uint8{1, 0, 0, 0, 0, 0})},
					{Coefficient: IntModP{1}, Exponents: NewBitPackedExponents([6]uint8{0, 1, 0, 0, 0, 0})},
					{Coefficient: IntModP{1}, Exponents: NewBitPackedExponents([6]uint8{0, 0, 1, 0, 0, 0})},
					{Coefficient: IntModP{1}, Exponents: NewBitPackedExponents([6]uint8{0, 0, 0, 1, 0, 0})},
					{Coefficient: IntModP{1}, Exponents: NewBitPackedExponents([6]uint8{0, 0, 0, 0, 1, 0})},
					{Coefficient: IntModP{1}, Exponents: NewBitPackedExponents([6]uint8{0, 0, 0, 0, 0, 1})},
				})
				p2 := NewPolynomial([]Term[IntModP, BitPackedExponents]{
					{Coefficient: IntModP{1}, Exponents: NewBitPackedExponents([6]uint8{1, 1, 0, 0, 0, 0})},
					{Coefficient: IntModP{1}, Exponents: NewBitPackedExponents([6]uint8{0, 1, 1, 0, 0, 0})},
					{Coefficient: IntModP{1}, Exponents: NewBitPackedExponents([6]uint8{0, 0, 1, 1, 0, 0})},
					{Coefficient: IntModP{1}, Exponents: NewBitPackedExponents([6]uint8{0, 0, 0, 1, 1, 0})},
					{Coefficient: IntModP{1}, Exponents: NewBitPackedExponents([6]uint8{0, 0, 0, 0, 1, 1})},
					{Coefficient: IntModP{1}, Exponents: NewBitPackedExponents([6]uint8{1, 0, 0, 0, 0, 1})},
				})
				p3 := NewPolynomial([]Term[IntModP, BitPackedExponents]{
					{Coefficient: IntModP{1}, Exponents: NewBitPackedExponents([6]uint8{1, 1, 1, 0, 0, 0})},
					{Coefficient: IntModP{1}, Exponents: NewBitPackedExponents([6]uint8{0, 1, 1, 1, 0, 0})},
					{Coefficient: IntModP{1}, Exponents: NewBitPackedExponents([6]uint8{0, 0, 1, 1, 1, 0})},
					{Coefficient: IntModP{1}, Exponents: NewBitPackedExponents([6]uint8{0, 0, 0, 1, 1, 1})},
					{Coefficient: IntModP{1}, Exponents: NewBitPackedExponents([6]uint8{1, 0, 0, 0, 1, 1})},
					{Coefficient: IntModP{1}, Exponents: NewBitPackedExponents([6]uint8{1, 1, 0, 0, 0, 1})},
				})
				p4 := NewPolynomial([]Term[IntModP, BitPackedExponents]{
					{Coefficient: IntModP{1}, Exponents: NewBitPackedExponents([6]uint8{1, 1, 1, 1, 0, 0})},
					{Coefficient: IntModP{1}, Exponents: NewBitPackedExponents([6]uint8{0, 1, 1, 1, 1, 0})},
					{Coefficient: IntModP{1}, Exponents: NewBitPackedExponents([6]uint8{0, 0, 1, 1, 1, 1})},
					{Coefficient: IntModP{1}, Exponents: NewBitPackedExponents([6]uint8{1, 0, 0, 1, 1, 1})},
					{Coefficient: IntModP{1}, Exponents: NewBitPackedExponents([6]uint8{1, 1, 0, 0, 1, 1})},
					{Coefficient: IntModP{1}, Exponents: NewBitPackedExponents([6]uint8{1, 1, 1, 0, 0, 1})},
				})
				p5 := NewPolynomial([]Term[IntModP, BitPackedExponents]{
					{Coefficient: IntModP{1}, Exponents: NewBitPackedExponents([6]uint8{1, 1, 1, 1, 1, 0})},
					{Coefficient: IntModP{1}, Exponents: NewBitPackedExponents([6]uint8{0, 1, 1, 1, 1, 1})},
					{Coefficient: IntModP{1}, Exponents: NewBitPackedExponents([6]uint8{1, 0, 1, 1, 1, 1})},
					{Coefficient: IntModP{1}, Exponents: NewBitPackedExponents([6]uint8{1, 1, 0, 1, 1, 1})},
					{Coefficient: IntModP{1}, Exponents: NewBitPackedExponents([6]uint8{1, 1, 1, 0, 1, 1})},
					{Coefficient: IntModP{1}, Exponents: NewBitPackedExponents([6]uint8{1, 1, 1, 1, 0, 1})},
				})
				p6 := NewPolynomial([]Term[IntModP, BitPackedExponents]{
					{Coefficient: IntModP{1}, Exponents: NewBitPackedExponents([6]uint8{1, 1, 1, 1, 1, 1})},
					{Coefficient: IntModP{modulus - 1}, Exponents: NewBitPackedExponents([6]uint8{0, 0, 0, 0, 0, 0})},
				})
				start := []Polynomial[IntModP, BitPackedExponents]{p1, p2, p3, p4, p5, p6}
				for i := 0; i < 10; i++ {
					basis := NaiveGrobnerBasis(start)
					fmt.Printf("Iteration %d: complete\n", i)
					if i == 9 {
						fmt.Println("Final Grobner Basis:")
						for _, poly := range basis {
							fmt.Printf("%v\n\n", poly)
						}
					}
				}
			}
		}
	}
}
