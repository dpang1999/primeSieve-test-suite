package finiteGrobner

import (
	"algos/helpers"
	"fmt"
	"hash/fnv"
	"slices"
	"sort"
)

type TermOrder int

const (
	Lex TermOrder = iota
	GrLex
	RevLex
)

var termOrder TermOrder = Lex
var modulus uint32 = 13

type Term struct {
	Coefficient uint32
	Exponents   []int
}

type Polynomial struct {
	Terms []Term
}

func (t Term) Compare(other Term) int {
	switch termOrder {
	case Lex:
		for i := range t.Exponents {
			if t.Exponents[i] != other.Exponents[i] {
				return t.Exponents[i] - other.Exponents[i]
			}
		}
		return 0
	case GrLex:
		tDeg := 0
		oDeg := 0
		for i := range t.Exponents {
			tDeg += t.Exponents[i]
			oDeg += other.Exponents[i]
		}
		if tDeg != oDeg {
			return tDeg - oDeg
		}
		for i := range t.Exponents {
			if t.Exponents[i] != other.Exponents[i] {
				return t.Exponents[i] - other.Exponents[i]
			}
		}
		return 0
	case RevLex:
		tDeg := 0
		oDeg := 0
		for i := range t.Exponents {
			tDeg += t.Exponents[i]
			oDeg += other.Exponents[i]
		}
		if tDeg != oDeg {
			return tDeg - oDeg
		}
		for i := len(t.Exponents) - 1; i >= 0; i-- {
			if t.Exponents[i] != other.Exponents[i] {
				return t.Exponents[i] - other.Exponents[i]
			}
		}
		return 0
	default:
		return 0
	}
}

func NewPolynomial(terms []Term) Polynomial {
	filtered := make([]Term, 0, len(terms))
	for _, t := range terms {
		if t.Coefficient != 0 {
			filtered = append(filtered, t)
		}
	}
	sort.Slice(filtered, func(i, j int) bool {
		return filtered[i].Compare(filtered[j]) > 0
	})
	return Polynomial{Terms: filtered}
}

func (p Polynomial) Add(other Polynomial) Polynomial {
	result := append([]Term{}, p.Terms...)
	for _, t := range other.Terms {
		found := false
		for i := range result {
			if slices.Equal(result[i].Exponents, t.Exponents) {
				result[i].Coefficient = (result[i].Coefficient + t.Coefficient) % modulus
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

func (p Polynomial) Subtract(other Polynomial) Polynomial {
	result := append([]Term{}, p.Terms...)
	for _, t := range other.Terms {
		found := false
		for i := range result {
			if slices.Equal(result[i].Exponents, t.Exponents) {
				result[i].Coefficient = (modulus + result[i].Coefficient - t.Coefficient) % modulus
				found = true
				break
			}
		}
		if !found {
			t2 := t
			t2.Coefficient = (modulus + 0 - t.Coefficient) % modulus
			result = append(result, t2)
		}
	}
	return NewPolynomial(result)
}

func (p Polynomial) MakeMonic() Polynomial {
	if len(p.Terms) == 0 {
		return p
	}

	leadCoeff := p.Terms[0].Coefficient
	if leadCoeff == 0 {
		return p
	}

	inv := modInverse(leadCoeff, modulus)
	terms := make([]Term, len(p.Terms))
	for i, t := range p.Terms {
		exps := append([]int(nil), t.Exponents...)
		terms[i] = Term{
			Coefficient: (t.Coefficient * inv) % modulus,
			Exponents:   exps,
		}
	}

	return NewPolynomial(terms)
}

func (p Polynomial) MultiplyByTerm(term Term) Polynomial {
	terms := make([]Term, len(p.Terms))
	for i, t := range p.Terms {
		newExps := make([]int, len(t.Exponents))
		for j := range t.Exponents {
			newExps[j] = t.Exponents[j] + term.Exponents[j]
		}
		terms[i] = Term{
			Coefficient: (t.Coefficient * term.Coefficient) % modulus,
			Exponents:   newExps,
		}
	}
	return NewPolynomial(terms)
}

func modInverse(a, m uint32) uint32 {
	var m0, x0, x1 = m, int64(0), int64(1)
	var aa, mm = int64(a), int64(m)
	if m == 1 {
		return 0
	}
	for aa > 1 {
		q := aa / mm
		t := mm
		mm = aa % mm
		aa = t
		tmp := x0
		x0 = x1 - q*x0
		x1 = tmp
	}
	if x1 < 0 {
		x1 += int64(m0)
	}
	return uint32(x1 % int64(m0))
}

// Fast hash for polynomial uniqueness
func polynomialHash(terms []Term) uint64 {
	h := fnv.New64a()
	for _, term := range terms {
		// Write coefficient
		var coeffBytes [4]byte
		coeff := term.Coefficient
		coeffBytes[0] = byte(coeff >> 24)
		coeffBytes[1] = byte(coeff >> 16)
		coeffBytes[2] = byte(coeff >> 8)
		coeffBytes[3] = byte(coeff)
		h.Write(coeffBytes[:])
		// Write exponents
		for _, exp := range term.Exponents {
			var expBytes [4]byte
			e := uint32(exp)
			expBytes[0] = byte(e >> 24)
			expBytes[1] = byte(e >> 16)
			expBytes[2] = byte(e >> 8)
			expBytes[3] = byte(e)
			h.Write(expBytes[:])
		}
	}
	return h.Sum64()
}

func (p Polynomial) Reduce(divisors []Polynomial) Polynomial {
	result := p
	remainder := []Term{}
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
			if canReduce(lead.Exponents, divLead.Exponents) {
				coeff := (lead.Coefficient * modInverse(divLead.Coefficient, modulus)) % modulus
				exps := make([]int, len(lead.Exponents))
				for i := range exps {
					exps[i] = lead.Exponents[i] - divLead.Exponents[i]
				}
				reductionTerm := Term{Coefficient: coeff, Exponents: exps}
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

func SPolynomial(p1 Polynomial, p2 Polynomial) Polynomial {
	lead1 := p1.Terms[0]
	lead2 := p2.Terms[0]
	lcmExps := make([]int, len(lead1.Exponents))
	for i := range lcmExps {
		lcmExps[i] = max(lead1.Exponents[i], lead2.Exponents[i])
	}
	scale1 := make([]int, len(lcmExps))
	scale2 := make([]int, len(lcmExps))
	for i := range lcmExps {
		scale1[i] = lcmExps[i] - lead1.Exponents[i]
		scale2[i] = lcmExps[i] - lead2.Exponents[i]
	}
	scaled1 := p1.MultiplyByTerm(Term{Coefficient: lead2.Coefficient, Exponents: scale1})
	scaled2 := p2.MultiplyByTerm(Term{Coefficient: lead1.Coefficient, Exponents: scale2})
	result := scaled1.Subtract(scaled2)
	return result
}

func NaiveGrobnerBasis(polys []Polynomial) []Polynomial {
	// ...existing code...
	basis := append([]Polynomial{}, polys...)
	basisSet := make(map[uint64]struct{})

	// Initialize basis set
	for _, poly := range basis {
		key := polynomialHash(poly.Terms)
		basisSet[key] = struct{}{}
	}

	// Create initial pairs
	type pair struct {
		i, j int
	}
	var pairs []pair
	for i := 0; i < len(basis); i++ {
		for j := i + 1; j < len(basis); j++ {
			pairs = append(pairs, pair{i, j})
		}
	}

	// Process pairs until none remain
	for len(pairs) > 0 {
		// Remove first pair
		p := pairs[0]
		pairs = pairs[1:]

		sPoly := SPolynomial(basis[p.i], basis[p.j])
		reduced := sPoly.Reduce(basis)
		key := polynomialHash(reduced.Terms)

		if len(reduced.Terms) > 0 {
			if _, ok := basisSet[key]; !ok {
				basisSet[key] = struct{}{}
				newIdx := len(basis)
				basis = append(basis, reduced)

				// Add new pairs with the new polynomial
				for k := 0; k < newIdx; k++ {
					pairs = append(pairs, pair{k, newIdx})
				}
			}
		}
	}
	//fmt.Printf("Basis before self-reduction: %v\n", basis)
	// Reduce basis by self
	reducedBasis := []Polynomial{}
	for _, poly := range basis {
		basisExcl := make([]Polynomial, 0, len(basis)-1)
		for _, p := range basis {
			if !polynomialsEqual(p, poly) {
				basisExcl = append(basisExcl, p)
			}
		}
		// compare length of basisExcl to basis
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

// Check if two exponent slices are equal
func equalExponents(a, b []int) bool {
	if slices.Equal(a, b) {
		return true
	} else {
		return false
	}

	/* if len(a) != len(b) {
		return false
	}
	for i := range a {
		if a[i] != b[i] {
			return false
		}
	}
	return true
	*/
}

func canReduce(a, b []int) bool {
	for i := range a {
		if a[i] < b[i] {
			return false
		}
	}
	return true
}

func polynomialsEqual(a, b Polynomial) bool {
	if len(a.Terms) != len(b.Terms) {
		return false
	}
	for i := range a.Terms {
		if a.Terms[i].Coefficient != b.Terms[i].Coefficient || !equalExponents(a.Terms[i].Exponents, b.Terms[i].Exponents) {
			return false
		}
	}
	return true
}

func TestFiniteGrobner(polyNum int, orderInt int, mod uint32) {
	// let mode == 0 be for testing
	mode := 0
	if mode != 0 {
		rand := helpers.NewLCG(12345, 1345, 16645, 1013904)
		if polyNum <= 0 {
			polyNum = 3
		}
		if orderInt < 0 || orderInt > 2 {
			orderInt = 0
		}
		if mod == 0 {
			mod = 13
		}
		modulus = mod
		numTerms := 3
		numExponents := 3
		termOrder = Lex
		var polys []Polynomial
		for p := 0; p < polyNum; p++ {
			var terms []Term
			for t := 0; t < numTerms; t++ {
				coeff := uint32(rand.NextInt()) % modulus
				exps := make([]int, numExponents)
				for i := range exps {
					exps[i] = rand.NextInt() % 4
				}
				terms = append(terms, Term{Coefficient: coeff, Exponents: exps})
			}
			polys = append(polys, NewPolynomial(terms))
		}
		basis := NaiveGrobnerBasis(polys)
		fmt.Println("Grobner Basis:")
		for i, poly := range basis {
			fmt.Printf("G%d: ", i)
			for _, term := range poly.Terms {
				fmt.Printf("%d*%v ", term.Coefficient, term.Exponents)
			}
			fmt.Println()
		}
	} else {
		n := polyNum
		if n <= 0 {
			n = 4
		}
		modulus = 7
		termOrder = Lex

		if n == 4 {
			// Cyclic 4
			fmt.Println("Go specialized finite coeff vec exponent cyclic 4")
			// a + b + c + d
			q1 := NewPolynomial([]Term{
				{Coefficient: 1, Exponents: []int{1, 0, 0, 0}},
				{Coefficient: 1, Exponents: []int{0, 1, 0, 0}},
				{Coefficient: 1, Exponents: []int{0, 0, 1, 0}},
				{Coefficient: 1, Exponents: []int{0, 0, 0, 1}},
			})
			// ab + bc + cd + ad
			q2 := NewPolynomial([]Term{
				{Coefficient: 1, Exponents: []int{1, 1, 0, 0}},
				{Coefficient: 1, Exponents: []int{0, 1, 1, 0}},
				{Coefficient: 1, Exponents: []int{0, 0, 1, 1}},
				{Coefficient: 1, Exponents: []int{1, 0, 0, 1}},
			})
			// abc + bcd + cda + dab
			q3 := NewPolynomial([]Term{
				{Coefficient: 1, Exponents: []int{1, 1, 1, 0}},
				{Coefficient: 1, Exponents: []int{0, 1, 1, 1}},
				{Coefficient: 1, Exponents: []int{1, 0, 1, 1}},
				{Coefficient: 1, Exponents: []int{1, 1, 0, 1}},
			})
			// abcd - 1
			q4 := NewPolynomial([]Term{
				{Coefficient: 1, Exponents: []int{1, 1, 1, 1}},
				{Coefficient: modulus - 1, Exponents: []int{0, 0, 0, 0}},
			})
			start := []Polynomial{q1, q2, q3, q4}
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
		} else if n == 5 {
			// Cyclic 5
			fmt.Println("Go specialized finite coeff vec exponent cyclic 5")
			// f1 = x0 + x1 + x2 + x3 + x4
			p1 := NewPolynomial([]Term{
				{Coefficient: 1, Exponents: []int{1, 0, 0, 0, 0}},
				{Coefficient: 1, Exponents: []int{0, 1, 0, 0, 0}},
				{Coefficient: 1, Exponents: []int{0, 0, 1, 0, 0}},
				{Coefficient: 1, Exponents: []int{0, 0, 0, 1, 0}},
				{Coefficient: 1, Exponents: []int{0, 0, 0, 0, 1}},
			})
			// f2 = x0x1 + x1x2 + x2x3 + x3x4 + x4x0
			p2 := NewPolynomial([]Term{
				{Coefficient: 1, Exponents: []int{1, 1, 0, 0, 0}},
				{Coefficient: 1, Exponents: []int{0, 1, 1, 0, 0}},
				{Coefficient: 1, Exponents: []int{0, 0, 1, 1, 0}},
				{Coefficient: 1, Exponents: []int{0, 0, 0, 1, 1}},
				{Coefficient: 1, Exponents: []int{1, 0, 0, 0, 1}},
			})
			// f3 = x0x1x2 + x1x2x3 + x2x3x4 + x3x4x0 + x4x0x1
			p3 := NewPolynomial([]Term{
				{Coefficient: 1, Exponents: []int{1, 1, 1, 0, 0}},
				{Coefficient: 1, Exponents: []int{0, 1, 1, 1, 0}},
				{Coefficient: 1, Exponents: []int{0, 0, 1, 1, 1}},
				{Coefficient: 1, Exponents: []int{1, 0, 0, 1, 1}},
				{Coefficient: 1, Exponents: []int{1, 1, 0, 0, 1}},
			})
			// f4 = x0x1x2x3 + x1x2x3x4 + x2x3x4x0 + x3x4x0x1 + x4x0x1x2
			p4 := NewPolynomial([]Term{
				{Coefficient: 1, Exponents: []int{1, 1, 1, 1, 0}},
				{Coefficient: 1, Exponents: []int{0, 1, 1, 1, 1}},
				{Coefficient: 1, Exponents: []int{1, 0, 1, 1, 1}},
				{Coefficient: 1, Exponents: []int{1, 1, 0, 1, 1}},
				{Coefficient: 1, Exponents: []int{1, 1, 1, 0, 1}},
			})
			// f5 = x0*x1*x2*x3*x4 - 1
			p5 := NewPolynomial([]Term{
				{Coefficient: 1, Exponents: []int{1, 1, 1, 1, 1}},
				{Coefficient: modulus - 1, Exponents: []int{0, 0, 0, 0, 0}},
			})
			start := []Polynomial{p1, p2, p3, p4, p5}
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
		} else if n == 6 {
			// Cyclic 6
			fmt.Println("Go specialized finite coeff vec exponent cyclic 6")
			// f1 = x0 + x1 + x2 + x3 + x4 + x5
			p1 := NewPolynomial([]Term{
				{Coefficient: 1, Exponents: []int{1, 0, 0, 0, 0, 0}},
				{Coefficient: 1, Exponents: []int{0, 1, 0, 0, 0, 0}},
				{Coefficient: 1, Exponents: []int{0, 0, 1, 0, 0, 0}},
				{Coefficient: 1, Exponents: []int{0, 0, 0, 1, 0, 0}},
				{Coefficient: 1, Exponents: []int{0, 0, 0, 0, 1, 0}},
				{Coefficient: 1, Exponents: []int{0, 0, 0, 0, 0, 1}},
			})
			// f2 = x0x1 + x1x2 + x2x3 + x3x4 + x4x5 + x5x0
			p2 := NewPolynomial([]Term{
				{Coefficient: 1, Exponents: []int{1, 1, 0, 0, 0, 0}},
				{Coefficient: 1, Exponents: []int{0, 1, 1, 0, 0, 0}},
				{Coefficient: 1, Exponents: []int{0, 0, 1, 1, 0, 0}},
				{Coefficient: 1, Exponents: []int{0, 0, 0, 1, 1, 0}},
				{Coefficient: 1, Exponents: []int{0, 0, 0, 0, 1, 1}},
				{Coefficient: 1, Exponents: []int{1, 0, 0, 0, 0, 1}},
			})
			// f3 = x0x1x2 + x1x2x3 + x2x3x4 + x3x4x5 + x4x5x0 + x5x0x1
			p3 := NewPolynomial([]Term{
				{Coefficient: 1, Exponents: []int{1, 1, 1, 0, 0, 0}},
				{Coefficient: 1, Exponents: []int{0, 1, 1, 1, 0, 0}},
				{Coefficient: 1, Exponents: []int{0, 0, 1, 1, 1, 0}},
				{Coefficient: 1, Exponents: []int{0, 0, 0, 1, 1, 1}},
				{Coefficient: 1, Exponents: []int{1, 0, 0, 0, 1, 1}},
				{Coefficient: 1, Exponents: []int{1, 1, 0, 0, 0, 1}},
			})
			// f4 = x0x1x2x3 + x1x2x3x4 + x2x3x4x5 + x3x4x5x0 + x4x5x0x1 + x5x0x1x2
			p4 := NewPolynomial([]Term{
				{Coefficient: 1, Exponents: []int{1, 1, 1, 1, 0, 0}},
				{Coefficient: 1, Exponents: []int{0, 1, 1, 1, 1, 0}},
				{Coefficient: 1, Exponents: []int{0, 0, 1, 1, 1, 1}},
				{Coefficient: 1, Exponents: []int{1, 0, 0, 1, 1, 1}},
				{Coefficient: 1, Exponents: []int{1, 1, 0, 0, 1, 1}},
				{Coefficient: 1, Exponents: []int{1, 1, 1, 0, 0, 1}},
			})
			// f5 = x0x1x2x3x4 + x1x2x3x4x5 + x2x3x4x5x0 + x3x4x5x0x1 + x4x5x0x1x2 + x5x0x1x2x3
			p5 := NewPolynomial([]Term{
				{Coefficient: 1, Exponents: []int{1, 1, 1, 1, 1, 0}},
				{Coefficient: 1, Exponents: []int{0, 1, 1, 1, 1, 1}},
				{Coefficient: 1, Exponents: []int{1, 0, 1, 1, 1, 1}},
				{Coefficient: 1, Exponents: []int{1, 1, 0, 1, 1, 1}},
				{Coefficient: 1, Exponents: []int{1, 1, 1, 0, 1, 1}},
				{Coefficient: 1, Exponents: []int{1, 1, 1, 1, 0, 1}},
			})
			// f6 = x0x1x2x3x4x5 - 1
			p6 := NewPolynomial([]Term{
				{Coefficient: 1, Exponents: []int{1, 1, 1, 1, 1, 1}},
				{Coefficient: modulus - 1, Exponents: []int{0, 0, 0, 0, 0, 0}},
			})
			start := []Polynomial{p1, p2, p3, p4, p5, p6}
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
			fmt.Println("Invalid test number")
		}
	}
}
