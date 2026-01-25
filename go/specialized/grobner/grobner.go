package grobner

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
	Exponents   []int
}

type Polynomial struct {
	Terms []Term
}

func (t Term) Compare(other Term, order TermOrder) int {
	switch order {
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

func NewPolynomial(terms []Term, order TermOrder) Polynomial {
	// Remove near-zero coefficients and round
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
			if equalExponents(result[i].Exponents, t.Exponents) {
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
			if equalExponents(result[i].Exponents, t.Exponents) {
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
		newExps := make([]int, len(t.Exponents))
		for j := range t.Exponents {
			newExps[j] = t.Exponents[j] + term.Exponents[j]
		}
		terms[i] = Term{
			Coefficient: t.Coefficient * term.Coefficient,
			Exponents:   newExps,
		}
	}
	return NewPolynomial(terms, order)
}

func (p Polynomial) Reduce(divisors []Polynomial, order TermOrder) Polynomial {
	//fmt.Printf("Reducing polynomial %v by divisors %v\n", p.Terms, divisors)
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
			if canReduce(lead.Exponents, divLead.Exponents) {
				coeff := lead.Coefficient / divLead.Coefficient
				exps := make([]int, len(lead.Exponents))
				for i := range exps {
					exps[i] = lead.Exponents[i] - divLead.Exponents[i]
				}
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
	//fmt.Printf("Final reduced polynomial: %v\n", result.Terms)
	return NewPolynomial(result.Terms, order)
}

func SPolynomial(p1 Polynomial, p2 Polynomial, order TermOrder) Polynomial {
	//fmt.Printf("Computing S-polynomial for %v and %v\n", p1.Terms, p2.Terms)
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
	scaled1 := p1.MultiplyByTerm(Term{Coefficient: 1.0, Exponents: scale1}, order)
	scaled2 := p2.MultiplyByTerm(Term{Coefficient: 1.0, Exponents: scale2}, order)
	result := scaled1.Subtract(scaled2, order)
	return result
}

func NaiveGrobnerBasis(polys []Polynomial, order TermOrder) []Polynomial {
	//fmt.Printf("Starting Grobner basis computation for polynomials: %v\n", polys)
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
		// print basis
		//fmt.Printf("Current basis: %v\n", basis)

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

// Check if two exponent slices are equal
func equalExponents(a, b []int) bool {
	if len(a) != len(b) {
		return false
	}
	for i := range a {
		if a[i] != b[i] {
			return false
		}
	}
	return true
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

func TestGrobner(polyNum int, orderInt int) {
	// let mode == 0 be for testing
	mode := 1
	if mode != 0 {
		// arg1 = # of polynomials
		// arg2 = term order (0=Lex,1=GrLex,2=RevLex)
		rand := helpers.NewLCG(12345, 1345, 16645, 1013904)
		if polyNum <= 0 {
			polyNum = 3 // default number of polynomials
		}
		if orderInt < 0 || orderInt > 2 {
			orderInt = 0 // default to Lex
		}
		numTerms := 3 // fixed number of terms per polynomial
		numExponents := 3
		order := TermOrder(orderInt)
		var polys []Polynomial
		for p := 0; p < polyNum; p++ {
			var terms []Term
			for t := 0; t < numTerms; t++ {
				coeff := (rand.NextDouble()*2 - 1)
				exps := make([]int, numExponents)
				for i := range exps {
					exps[i] = rand.NextInt() % 4 // exponents 0 to 3
				}
				terms = append(terms, Term{Coefficient: coeff, Exponents: exps})
			}
			polys = append(polys, NewPolynomial(terms, order))
		}
		basis := NaiveGrobnerBasis(polys, order)
		fmt.Println("Grobner Basis:")
		for i, poly := range basis {
			fmt.Printf("G%d: ", i)
			for _, term := range poly.Terms {
				fmt.Printf("%.2f*%v ", term.Coefficient, term.Exponents)
			}
			fmt.Println()
		}
	} else {

		fmt.Println("--- Grobner Basis Test ---")

		// x^3 + y^3 + z^3
		p1 := NewPolynomial([]Term{
			{Coefficient: 1, Exponents: []int{3, 0, 0}},
			{Coefficient: 1, Exponents: []int{0, 3, 0}},
			{Coefficient: 1, Exponents: []int{0, 0, 3}},
		}, Lex)

		// xy + yz + xz
		p2 := NewPolynomial([]Term{
			{Coefficient: 1, Exponents: []int{1, 0, 1}},
			{Coefficient: 1, Exponents: []int{0, 1, 1}},
			{Coefficient: 1, Exponents: []int{1, 1, 0}},
		}, Lex)

		//x+y+z
		p3 := NewPolynomial([]Term{
			{Coefficient: 1, Exponents: []int{1, 0, 0}},
			{Coefficient: 1, Exponents: []int{0, 1, 0}},
			{Coefficient: 1, Exponents: []int{0, 0, 1}},
		}, Lex)
		fmt.Printf("Input polynomials:\nP1: %v\nP2: %v\nP3: %v\n", p1.Terms, p2.Terms, p3.Terms)
		basis := NaiveGrobnerBasis([]Polynomial{p1, p2, p3}, Lex)

		fmt.Println("Grobner Basis:")
		for i, poly := range basis {
			fmt.Printf("G%d: ", i)
			for _, term := range poly.Terms {
				fmt.Printf("%.2f*%v ", term.Coefficient, term.Exponents)
			}
			fmt.Println()
		}
		fmt.Println("--- End Grobner Basis Test ---")
	}
}
