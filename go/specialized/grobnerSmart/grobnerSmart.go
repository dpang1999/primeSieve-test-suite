package grobnerSmart

import (
	"algos/helpers"
	"fmt"
	"math"
	"sort"
)

type TermOrder int

const (
	Lex TermOrder = iota
	GrLex
	RevLex
)

type Term struct {
	Coefficient float64
	Exponents   uint64 // Bitpacked: [63..48]=degree, [47..40]=e0, [39..32]=e1, ... [7..0]=e5
}

type Polynomial struct {
	Terms []Term
}

func (t Term) Compare(other Term, order TermOrder) int {
	switch order {
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

func NewPolynomial(terms []Term, order TermOrder) Polynomial {
	filtered := make([]Term, 0, len(terms))
	for _, t := range terms {
		if math.Abs(t.Coefficient) > 1e-2 {
			t.Coefficient = math.Round(t.Coefficient*1e5) / 1e5
			filtered = append(filtered, t)
		}
	}
	sort.Slice(filtered, func(i, j int) bool {
		return filtered[i].Compare(filtered[j], order) > 0
	})
	return Polynomial{Terms: filtered}
}

func (p Polynomial) Add(other Polynomial, order TermOrder) Polynomial {
	result := append([]Term{}, p.Terms...)
	for _, t := range other.Terms {
		found := false
		for i := range result {
			if result[i].Exponents == t.Exponents {
				result[i].Coefficient += t.Coefficient
				found = true
				break
			}
		}
		if !found {
			result = append(result, t)
		}
	}
	return NewPolynomial(result, order)
}

func (p Polynomial) Subtract(other Polynomial, order TermOrder) Polynomial {
	result := append([]Term{}, p.Terms...)
	for _, t := range other.Terms {
		found := false
		for i := range result {
			if result[i].Exponents == t.Exponents {
				result[i].Coefficient -= t.Coefficient
				found = true
				break
			}
		}
		if !found {
			t2 := t
			t2.Coefficient = -t2.Coefficient
			result = append(result, t2)
		}
	}
	return NewPolynomial(result, order)
}

func (p Polynomial) MultiplyByTerm(term Term, order TermOrder) Polynomial {
	terms := make([]Term, len(p.Terms))
	for i, t := range p.Terms {
		terms[i] = Term{
			Coefficient: t.Coefficient * term.Coefficient,
			Exponents:   t.Exponents + term.Exponents,
		}
	}
	return NewPolynomial(terms, order)
}

func (p Polynomial) Reduce(divisors []Polynomial, order TermOrder) Polynomial {
	result := p
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
				coeff := lead.Coefficient / divLead.Coefficient
				exps := lead.Exponents - divLead.Exponents
				reductionTerm := Term{Coefficient: coeff, Exponents: exps}
				scaledDivisor := divisor.MultiplyByTerm(reductionTerm, order)
				result = result.Subtract(scaledDivisor, order)
				reduced = true
				break
			}
		}
		if !reduced {
			break
		}
	}
	return NewPolynomial(result.Terms, order)
}

func SPolynomial(p1 Polynomial, p2 Polynomial, order TermOrder) Polynomial {
	lead1 := p1.Terms[0]
	lead2 := p2.Terms[0]
	lcmExps := lcmPackedExponents(lead1.Exponents, lead2.Exponents)
	scale1 := lcmExps - lead1.Exponents
	scale2 := lcmExps - lead2.Exponents
	scaled1 := p1.MultiplyByTerm(Term{Coefficient: 1.0, Exponents: scale1}, order)
	scaled2 := p2.MultiplyByTerm(Term{Coefficient: 1.0, Exponents: scale2}, order)
	result := scaled1.Subtract(scaled2, order)
	return result
}

func NaiveGrobnerBasis(polys []Polynomial, order TermOrder) []Polynomial {
	basis := append([]Polynomial{}, polys...)
	basisSet := make(map[string]struct{})
	added := true
	for added {
		added = false
		n := len(basis)
		for i := 0; i < n; i++ {
			for j := i + 1; j < n; j++ {
				sPoly := SPolynomial(basis[i], basis[j], order)
				reduced := sPoly.Reduce(basis, order)
				key := fmt.Sprintf("%v", reduced.Terms)
				if len(reduced.Terms) > 0 {
					if _, ok := basisSet[key]; !ok {
						basisSet[key] = struct{}{}
						basis = append(basis, reduced)
						added = true
					}
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
		reduced := poly.Reduce(basisExcl, order)
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
	mode := 1
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
		order := TermOrder(orderInt)
		var polys []Polynomial
		for p := 0; p < polyNum; p++ {
			var terms []Term
			for t := 0; t < numTerms; t++ {
				coeff := rand.NextDouble()*2 - 1
				var exps [6]uint8
				for i := 0; i < numExponents; i++ {
					exps[i] = uint8(rand.NextInt() % 4)
				}
				terms = append(terms, Term{Coefficient: coeff, Exponents: packExponents(exps)})
			}
			polys = append(polys, NewPolynomial(terms, order))
		}
		basis := NaiveGrobnerBasis(polys, order)
		fmt.Println("Grobner Basis:")
		for i, poly := range basis {
			fmt.Printf("G%d: ", i)
			for _, term := range poly.Terms {
				fmt.Printf("%.2f*%v ", term.Coefficient, unpackExponents(term.Exponents))
			}
			fmt.Println()
		}
	} else {
		fmt.Println("--- Grobner Basis Test ---")
		// x^3 + y^3 + z^3
		p1 := NewPolynomial([]Term{
			{Coefficient: 1, Exponents: packExponents([6]uint8{3, 0, 0, 0, 0, 0})},
			{Coefficient: 1, Exponents: packExponents([6]uint8{0, 3, 0, 0, 0, 0})},
			{Coefficient: 1, Exponents: packExponents([6]uint8{0, 0, 3, 0, 0, 0})},
		}, Lex)
		// xy + yz + xz
		p2 := NewPolynomial([]Term{
			{Coefficient: 1, Exponents: packExponents([6]uint8{1, 0, 1, 0, 0, 0})},
			{Coefficient: 1, Exponents: packExponents([6]uint8{0, 1, 1, 0, 0, 0})},
			{Coefficient: 1, Exponents: packExponents([6]uint8{1, 1, 0, 0, 0, 0})},
		}, Lex)
		// x+y+z
		p3 := NewPolynomial([]Term{
			{Coefficient: 1, Exponents: packExponents([6]uint8{1, 0, 0, 0, 0, 0})},
			{Coefficient: 1, Exponents: packExponents([6]uint8{0, 1, 0, 0, 0, 0})},
			{Coefficient: 1, Exponents: packExponents([6]uint8{0, 0, 1, 0, 0, 0})},
		}, Lex)
		fmt.Printf("Input polynomials:\nP1: %v\nP2: %v\nP3: %v\n", p1.Terms, p2.Terms, p3.Terms)
		basis := NaiveGrobnerBasis([]Polynomial{p1, p2, p3}, Lex)
		fmt.Println("Grobner Basis:")
		for i, poly := range basis {
			fmt.Printf("G%d: ", i)
			for _, term := range poly.Terms {
				fmt.Printf("%.2f*%v ", term.Coefficient, unpackExponents(term.Exponents))
			}
			fmt.Println()
		}
		fmt.Println("--- End Grobner Basis Test ---")
	}
}
