package grobnerSmart

import (
	"algos/helpers"
	"fmt"
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
	Exponents   uint64 // Bitpacked: [63..48]=degree, [47..40]=e0, [39..32]=e1, ... [7..0]=e5
}

type Polynomial struct {
	Terms []Term
}

func (t Term) Compare(other Term) int {
	switch termOrder {
	case Lex:
		a := t.Exponents & 0x0000FFFFFFFFFFFF
		b := other.Exponents & 0x0000FFFFFFFFFFFF
		if a > b {
			return 1
		} else if a < b {
			return -1
		}
		return 0
	case GrLex:
		tDeg := int((t.Exponents >> 48) & 0xFFFF)
		oDeg := int((other.Exponents >> 48) & 0xFFFF)
		if tDeg != oDeg {
			return tDeg - oDeg
		}
		a := t.Exponents & 0x0000FFFFFFFFFFFF
		b := other.Exponents & 0x0000FFFFFFFFFFFF
		if a > b {
			return 1
		} else if a < b {
			return -1
		}
		return 0
	case RevLex:
		tDeg := int((t.Exponents >> 48) & 0xFFFF)
		oDeg := int((other.Exponents >> 48) & 0xFFFF)
		if tDeg != oDeg {
			return tDeg - oDeg
		}
		a := t.Exponents & 0x0000FFFFFFFFFFFF
		b := other.Exponents & 0x0000FFFFFFFFFFFF
		if a < b {
			return 1
		} else if a > b {
			return -1
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
			if result[i].Exponents == t.Exponents {
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
			if result[i].Exponents == t.Exponents {
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
		terms[i] = Term{
			Coefficient: (t.Coefficient * inv) % modulus,
			Exponents:   t.Exponents,
		}
	}

	return NewPolynomial(terms)
}

func (p Polynomial) MultiplyByTerm(term Term) Polynomial {
	terms := make([]Term, len(p.Terms))
	for i, t := range p.Terms {
		terms[i] = Term{
			Coefficient: t.Coefficient * term.Coefficient % modulus,
			Exponents:   t.Exponents + term.Exponents,
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
			if canReducePacked(lead.Exponents, divLead.Exponents) {
				coeff := (lead.Coefficient * modInverse(divLead.Coefficient, modulus)) % modulus
				exps := lead.Exponents - divLead.Exponents
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
	lcmExps := lcmPackedExponents(lead1.Exponents, lead2.Exponents)
	scale1 := lcmExps - lead1.Exponents
	scale2 := lcmExps - lead2.Exponents
	scaled1 := p1.MultiplyByTerm(Term{Coefficient: 1.0, Exponents: scale1})
	scaled2 := p2.MultiplyByTerm(Term{Coefficient: 1.0, Exponents: scale2})
	result := scaled1.Subtract(scaled2)
	return result
}

func NaiveGrobnerBasis(polys []Polynomial) []Polynomial {
	basis := append([]Polynomial{}, polys...)
	basisSet := make(map[string]struct{})

	for _, poly := range basis {
		key := fmt.Sprintf("%v", poly.Terms)
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
		sPoly := SPolynomial(basis[p.i], basis[p.j])
		reduced := sPoly.Reduce(basis)
		key := fmt.Sprintf("%v", reduced.Terms)
		if len(reduced.Terms) > 0 {
			if _, ok := basisSet[key]; !ok {
				basisSet[key] = struct{}{}
				newIdx := len(basis)
				basis = append(basis, reduced)
				for i := 0; i < newIdx; i++ {
					pairs = append(pairs, pair{i, newIdx})
				}
			}
		}
	}

	reducedBasis := []Polynomial{}
	for _, poly := range basis {
		basisExcl := make([]Polynomial, 0, len(basis)-1)
		for _, p := range basis {
			if !polynomialsEqual(p, poly) {
				basisExcl = append(basisExcl, p)
			}
		}
		reduced := poly.Reduce(basisExcl).MakeMonic()
		key := fmt.Sprintf("%v", reduced.Terms)
		if len(reduced.Terms) > 0 {
			found := false
			for _, rb := range reducedBasis {
				if fmt.Sprintf("%v", rb.Terms) == key {
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

// Bitpacked exponent helpers (6 variables, 8 bits each, top 16 bits for degree)
func packExponents(exps [6]uint8) uint64 {
	var packed uint64 = 0
	for i, e := range exps {
		shift := 40 - 8*i
		packed |= uint64(e) << shift
	}
	var degree uint16 = 0
	for _, e := range exps {
		degree += uint16(e)
	}
	packed |= uint64(degree) << 48
	return packed
}

func unpackExponents(packed uint64) [6]uint8 {
	var exps [6]uint8
	for i := 0; i < 6; i++ {
		shift := 40 - 8*i
		exps[i] = uint8((packed >> shift) & 0xFF)
	}
	return exps
}

// Compute the LCM of two packed exponent vectors directly, without unpacking
func lcmPackedExponents(a, b uint64) uint64 {
	var lcm uint64 = 0
	var degree uint16 = 0
	for i := 0; i < 6; i++ {
		shift := 40 - 8*i
		ea := (a >> shift) & 0xFF
		eb := (b >> shift) & 0xFF
		l := ea
		if eb > ea {
			l = eb
		}
		lcm |= l << shift
		degree += uint16(l)
	}
	lcm |= uint64(degree) << 48
	return lcm
}

func canReducePacked(a, b uint64) bool {
	for i := 0; i < 6; i++ {
		shift := 40 - 8*i
		ea := (a >> shift) & 0xFF
		eb := (b >> shift) & 0xFF
		if ea < eb {
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
		if a.Terms[i].Coefficient != b.Terms[i].Coefficient || a.Terms[i].Exponents != b.Terms[i].Exponents {
			return false
		}
	}
	return true
}

func TestGrobnerSmart(polyNum int, orderInt int) {
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
		numTerms := 3
		numExponents := 3
		var polys []Polynomial
		for p := 0; p < polyNum; p++ {
			var terms []Term
			for t := 0; t < numTerms; t++ {
				coeff := uint32(rand.NextInt()) % modulus
				var exps [6]uint8
				for i := 0; i < numExponents; i++ {
					exps[i] = uint8(rand.NextInt() % 4)
				}
				terms = append(terms, Term{Coefficient: coeff, Exponents: packExponents(exps)})
			}
			polys = append(polys, NewPolynomial(terms))
		}
		basis := NaiveGrobnerBasis(polys)
		fmt.Println("Grobner Basis:")
		for i, poly := range basis {
			fmt.Printf("G%d: ", i)
			for _, term := range poly.Terms {
				fmt.Printf("%d*%v ", term.Coefficient, unpackExponents(term.Exponents))
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
			fmt.Println("Go specialized finite coeff bitpacked exponent cyclic 4")
			// a + b + c + d
			q1 := NewPolynomial([]Term{
				{Coefficient: 1, Exponents: packExponents([6]uint8{1, 0, 0, 0, 0, 0})},
				{Coefficient: 1, Exponents: packExponents([6]uint8{0, 1, 0, 0, 0, 0})},
				{Coefficient: 1, Exponents: packExponents([6]uint8{0, 0, 1, 0, 0, 0})},
				{Coefficient: 1, Exponents: packExponents([6]uint8{0, 0, 0, 1, 0, 0})},
			})
			// ab + bc + cd + ad
			q2 := NewPolynomial([]Term{
				{Coefficient: 1, Exponents: packExponents([6]uint8{1, 1, 0, 0, 0, 0})},
				{Coefficient: 1, Exponents: packExponents([6]uint8{0, 1, 1, 0, 0, 0})},
				{Coefficient: 1, Exponents: packExponents([6]uint8{0, 0, 1, 1, 0, 0})},
				{Coefficient: 1, Exponents: packExponents([6]uint8{1, 0, 0, 1, 0, 0})},
			})
			// abc + bcd + cda + dab
			q3 := NewPolynomial([]Term{
				{Coefficient: 1, Exponents: packExponents([6]uint8{1, 1, 1, 0, 0, 0})},
				{Coefficient: 1, Exponents: packExponents([6]uint8{0, 1, 1, 1, 0, 0})},
				{Coefficient: 1, Exponents: packExponents([6]uint8{1, 0, 1, 1, 0, 0})},
				{Coefficient: 1, Exponents: packExponents([6]uint8{1, 1, 0, 1, 0, 0})},
			})
			// abcd - 1
			q4 := NewPolynomial([]Term{
				{Coefficient: 1, Exponents: packExponents([6]uint8{1, 1, 1, 1, 0, 0})},
				{Coefficient: modulus - 1, Exponents: packExponents([6]uint8{0, 0, 0, 0, 0, 0})},
			})
			start := []Polynomial{q1, q2, q3, q4}
			for i := 0; i < 10; i++ {
				basis := NaiveGrobnerBasis(start)
				fmt.Printf("Iteration %d: complete\n", i)
				if i == 9 {
					fmt.Println("Final Grobner Basis:")
					for i, poly := range basis {
						fmt.Printf("G%d: ", i)
						for _, term := range poly.Terms {
							fmt.Printf("%d*%v ", term.Coefficient, unpackExponents(term.Exponents))
						}
						fmt.Println()
					}
				}
			}
		} else if n == 5 {
			// Cyclic 5
			fmt.Println("Go specialized finite coeff bitpacked exponent cyclic 5")
			// f1 = x0 + x1 + x2 + x3 + x4
			p1 := NewPolynomial([]Term{
				{Coefficient: 1, Exponents: packExponents([6]uint8{1, 0, 0, 0, 0, 0})},
				{Coefficient: 1, Exponents: packExponents([6]uint8{0, 1, 0, 0, 0, 0})},
				{Coefficient: 1, Exponents: packExponents([6]uint8{0, 0, 1, 0, 0, 0})},
				{Coefficient: 1, Exponents: packExponents([6]uint8{0, 0, 0, 1, 0, 0})},
				{Coefficient: 1, Exponents: packExponents([6]uint8{0, 0, 0, 0, 1, 0})},
			})
			// f2 = x0x1 + x1x2 + x2x3 + x3x4 + x4x0
			p2 := NewPolynomial([]Term{
				{Coefficient: 1, Exponents: packExponents([6]uint8{1, 1, 0, 0, 0, 0})},
				{Coefficient: 1, Exponents: packExponents([6]uint8{0, 1, 1, 0, 0, 0})},
				{Coefficient: 1, Exponents: packExponents([6]uint8{0, 0, 1, 1, 0, 0})},
				{Coefficient: 1, Exponents: packExponents([6]uint8{0, 0, 0, 1, 1, 0})},
				{Coefficient: 1, Exponents: packExponents([6]uint8{1, 0, 0, 0, 1, 0})},
			})
			// f3 = x0x1x2 + x1x2x3 + x2x3x4 + x3x4x0 + x4x0x1
			p3 := NewPolynomial([]Term{
				{Coefficient: 1, Exponents: packExponents([6]uint8{1, 1, 1, 0, 0, 0})},
				{Coefficient: 1, Exponents: packExponents([6]uint8{0, 1, 1, 1, 0, 0})},
				{Coefficient: 1, Exponents: packExponents([6]uint8{0, 0, 1, 1, 1, 0})},
				{Coefficient: 1, Exponents: packExponents([6]uint8{1, 0, 0, 1, 1, 0})},
				{Coefficient: 1, Exponents: packExponents([6]uint8{1, 1, 0, 0, 1, 0})},
			})
			// f4 = x0x1x2x3 + x1x2x3x4 + x2x3x4x0 + x3x4x0x1 + x4x0x1x2
			p4 := NewPolynomial([]Term{
				{Coefficient: 1, Exponents: packExponents([6]uint8{1, 1, 1, 1, 0, 0})},
				{Coefficient: 1, Exponents: packExponents([6]uint8{0, 1, 1, 1, 1, 0})},
				{Coefficient: 1, Exponents: packExponents([6]uint8{1, 0, 1, 1, 1, 0})},
				{Coefficient: 1, Exponents: packExponents([6]uint8{1, 1, 0, 1, 1, 0})},
				{Coefficient: 1, Exponents: packExponents([6]uint8{1, 1, 1, 0, 1, 0})},
			})
			// f5 = x0*x1*x2*x3*x4 - 1
			p5 := NewPolynomial([]Term{
				{Coefficient: 1, Exponents: packExponents([6]uint8{1, 1, 1, 1, 1, 0})},
				{Coefficient: modulus - 1, Exponents: packExponents([6]uint8{0, 0, 0, 0, 0, 0})},
			})
			start := []Polynomial{p1, p2, p3, p4, p5}
			for i := 0; i < 10; i++ {
				basis := NaiveGrobnerBasis(start)
				fmt.Printf("Iteration %d: complete\n", i)
				if i == 9 {
					fmt.Println("Final Grobner Basis:")
					for i, poly := range basis {
						fmt.Printf("G%d: ", i)
						for _, term := range poly.Terms {
							fmt.Printf("%d*%v ", term.Coefficient, unpackExponents(term.Exponents))
						}
						fmt.Println()
					}
				}
			}
		} else if n == 6 {
			// Cyclic 6
			fmt.Println("Go specialized finite coeff bitpacked exponent cyclic 6")
			// f1 = x0 + x1 + x2 + x3 + x4 + x5
			p1 := NewPolynomial([]Term{
				{Coefficient: 1, Exponents: packExponents([6]uint8{1, 0, 0, 0, 0, 0})},
				{Coefficient: 1, Exponents: packExponents([6]uint8{0, 1, 0, 0, 0, 0})},
				{Coefficient: 1, Exponents: packExponents([6]uint8{0, 0, 1, 0, 0, 0})},
				{Coefficient: 1, Exponents: packExponents([6]uint8{0, 0, 0, 1, 0, 0})},
				{Coefficient: 1, Exponents: packExponents([6]uint8{0, 0, 0, 0, 1, 0})},
				{Coefficient: 1, Exponents: packExponents([6]uint8{0, 0, 0, 0, 0, 1})},
			})
			// f2 = x0x1 + x1x2 + x2x3 + x3x4 + x4x5 + x5x0
			p2 := NewPolynomial([]Term{
				{Coefficient: 1, Exponents: packExponents([6]uint8{1, 1, 0, 0, 0, 0})},
				{Coefficient: 1, Exponents: packExponents([6]uint8{0, 1, 1, 0, 0, 0})},
				{Coefficient: 1, Exponents: packExponents([6]uint8{0, 0, 1, 1, 0, 0})},
				{Coefficient: 1, Exponents: packExponents([6]uint8{0, 0, 0, 1, 1, 0})},
				{Coefficient: 1, Exponents: packExponents([6]uint8{0, 0, 0, 0, 1, 1})},
				{Coefficient: 1, Exponents: packExponents([6]uint8{1, 0, 0, 0, 0, 1})},
			})
			// f3 = x0x1x2 + x1x2x3 + x2x3x4 + x3x4x5 + x4x5x0 + x5x0x1
			p3 := NewPolynomial([]Term{
				{Coefficient: 1, Exponents: packExponents([6]uint8{1, 1, 1, 0, 0, 0})},
				{Coefficient: 1, Exponents: packExponents([6]uint8{0, 1, 1, 1, 0, 0})},
				{Coefficient: 1, Exponents: packExponents([6]uint8{0, 0, 1, 1, 1, 0})},
				{Coefficient: 1, Exponents: packExponents([6]uint8{0, 0, 0, 1, 1, 1})},
				{Coefficient: 1, Exponents: packExponents([6]uint8{1, 0, 0, 0, 1, 1})},
				{Coefficient: 1, Exponents: packExponents([6]uint8{1, 1, 0, 0, 0, 1})},
			})
			// f4 = x0x1x2x3 + x1x2x3x4 + x2x3x4x5 + x3x4x5x0 + x4x5x0x1 + x5x0x1x2
			p4 := NewPolynomial([]Term{
				{Coefficient: 1, Exponents: packExponents([6]uint8{1, 1, 1, 1, 0, 0})},
				{Coefficient: 1, Exponents: packExponents([6]uint8{0, 1, 1, 1, 1, 0})},
				{Coefficient: 1, Exponents: packExponents([6]uint8{0, 0, 1, 1, 1, 1})},
				{Coefficient: 1, Exponents: packExponents([6]uint8{1, 0, 0, 1, 1, 1})},
				{Coefficient: 1, Exponents: packExponents([6]uint8{1, 1, 0, 0, 1, 1})},
				{Coefficient: 1, Exponents: packExponents([6]uint8{1, 1, 1, 0, 0, 1})},
			})
			// f5 = x0x1x2x3x4 + x1x2x3x4x5 + x2x3x4x5x0 + x3x4x5x0x1 + x4x5x0x1x2 + x5x0x1x2x3
			p5 := NewPolynomial([]Term{
				{Coefficient: 1, Exponents: packExponents([6]uint8{1, 1, 1, 1, 1, 0})},
				{Coefficient: 1, Exponents: packExponents([6]uint8{0, 1, 1, 1, 1, 1})},
				{Coefficient: 1, Exponents: packExponents([6]uint8{1, 0, 1, 1, 1, 1})},
				{Coefficient: 1, Exponents: packExponents([6]uint8{1, 1, 0, 1, 1, 1})},
				{Coefficient: 1, Exponents: packExponents([6]uint8{1, 1, 1, 0, 1, 1})},
				{Coefficient: 1, Exponents: packExponents([6]uint8{1, 1, 1, 1, 0, 1})},
			})
			// f6 = x0x1x2x3x4x5 - 1
			p6 := NewPolynomial([]Term{
				{Coefficient: 1, Exponents: packExponents([6]uint8{1, 1, 1, 1, 1, 1})},
				{Coefficient: modulus - 1, Exponents: packExponents([6]uint8{0, 0, 0, 0, 0, 0})},
			})
			start := []Polynomial{p1, p2, p3, p4, p5, p6}
			for i := 0; i < 10; i++ {
				basis := NaiveGrobnerBasis(start)
				fmt.Printf("Iteration %d: complete\n", i)
				if i == 9 {
					fmt.Println("Final Grobner Basis:")
					for i, poly := range basis {
						fmt.Printf("G%d: ", i)
						for _, term := range poly.Terms {
							fmt.Printf("%d*%v ", term.Coefficient, unpackExponents(term.Exponents))
						}
						fmt.Println()
					}
				}
			}
		} else {
			fmt.Println("Invalid test number")
		}
	}
}
