package generic

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

func (t Term[N, E]) Compare(other Term[N, E], order TermOrder) int {
	switch order {
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
}](terms []Term[N, E], order TermOrder) Polynomial[N, E] {
	filtered := make([]Term[N, E], 0, len(terms))
	for _, t := range terms {
		if math.Abs(t.Coefficient.coerceToFloat()) > 1e-2 {
			t.Coefficient = t.Coefficient.coerceFromFloat(
				math.Round(t.Coefficient.coerceToFloat()*1e5) / 1e5)
			filtered = append(filtered, t)
		}
	}
	sort.Slice(filtered, func(i, j int) bool {
		return filtered[i].Compare(filtered[j], order) > 0
	})
	return Polynomial[N, E]{Terms: filtered}
}

func (p Polynomial[N, E]) Add(other Polynomial[N, E], order TermOrder) Polynomial[N, E] {
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
	return NewPolynomial(result, order)
}

func (p Polynomial[N, E]) Subtract(other Polynomial[N, E], order TermOrder) Polynomial[N, E] {
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
	return NewPolynomial(result, order)
}

func (p Polynomial[N, E]) MultiplyByTerm(term Term[N, E], order TermOrder) Polynomial[N, E] {
	terms := make([]Term[N, E], len(p.Terms))
	for i, t := range p.Terms {
		terms[i] = Term[N, E]{
			Coefficient: t.Coefficient.m(term.Coefficient),
			Exponents:   t.Exponents.add(term.Exponents),
		}
	}
	return NewPolynomial(terms, order)
}

func (p Polynomial[N, E]) Reduce(divisors []Polynomial[N, E], order TermOrder) Polynomial[N, E] {
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
			if lead.Exponents.canReduce(divLead.Exponents) {
				coeff := lead.Coefficient.d(divLead.Coefficient)
				exps := lead.Exponents.sub(divLead.Exponents)
				reductionTerm := Term[N, E]{Coefficient: coeff, Exponents: exps}
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

func SPolynomial[N interface {
	IField[N]
	IMath[N]
	ICopiable[N]
	IOrdered[N]
}, E interface {
	IExponents[E]
}](p1 Polynomial[N, E], p2 Polynomial[N, E], order TermOrder) Polynomial[N, E] {
	lead1 := p1.Terms[0]
	lead2 := p2.Terms[0]
	lcmExps := lead1.Exponents.lcm(lead2.Exponents)
	scale1 := lcmExps.sub(lead1.Exponents)
	scale2 := lcmExps.sub(lead2.Exponents)
	scaled1 := p1.MultiplyByTerm(Term[N, E]{Coefficient: lead1.Coefficient.one(), Exponents: scale1}, order)
	scaled2 := p2.MultiplyByTerm(Term[N, E]{Coefficient: lead2.Coefficient.one(), Exponents: scale2}, order)
	result := scaled1.Subtract(scaled2, order)
	return result
}

func NaiveGrobnerBasis[N interface {
	IField[N]
	IMath[N]
	ICopiable[N]
	IOrdered[N]
}, E interface {
	IExponents[E]
}](polys []Polynomial[N, E], order TermOrder) []Polynomial[N, E] {
	basis := append([]Polynomial[N, E]{}, polys...)
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
	reducedBasis := []Polynomial[N, E]{}
	for _, poly := range basis {
		basisExcl := make([]Polynomial[N, E], 0, len(basis)-1)
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
		rand := helpers.NewLCG(12345, 1345, 65, 17)
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
func TestGenGrobner(numPoly, numTerms, coeffType, expType, orderInt, modulus int) {
	fmt.Printf("Generating random basis: numPoly=%d, numTerms=%d, coeffType=%d, expType=%d, order=%d, modulus=%d\n",
		numPoly, numTerms, coeffType, expType, orderInt, modulus)
	if numPoly <= 0 {
		numPoly = 3
	}
	if numTerms <= 0 {
		numTerms = 3
	}
	if coeffType < 0 || coeffType > 2 {
		coeffType = 0
		// 0 = doublefield, 1 = intmodp, 2 = singlefield
	}
	if expType < 0 || expType > 1 {
		expType = 0
		// 0 = vecexponents, 1 = bitpackedexponents
	}
	if orderInt < 0 || orderInt > 2 {
		orderInt = 0
		// 0 = lex, 1 = grlex, 2 = revlex
	}
	if modulus <= 0 {
		modulus = 13
	}
	order := TermOrder(orderInt)
	rand := helpers.NewLCG(12345, 1345, 65, 17)

	switch {
	case coeffType == 0 && expType == 0:
		// DoubleField + VecExponents
		var polys []Polynomial[DoubleField, VecExponents]
		for i := 0; i < numPoly; i++ {
			var terms []Term[DoubleField, VecExponents]
			for j := 0; j < numTerms; j++ {
				coeff := DoubleField{Value: rand.NextDouble()}
				exps := NewVecExponents([]uint32{uint32(rand.NextInt() % 4), uint32(rand.NextInt() % 4), uint32(rand.NextInt() % 4)})
				terms = append(terms, Term[DoubleField, VecExponents]{Coefficient: coeff, Exponents: exps})
			}
			polys = append(polys, NewPolynomial(terms, order))
		}
		basis := NaiveGrobnerBasis(polys, order)
		fmt.Println("DoubleField + VecExponents:")
		for i, poly := range basis {
			fmt.Printf("G%d: %v\n", i, poly)
		}
	case coeffType == 0 && expType == 1:
		// DoubleField + BitPackedExponents
		var polys []Polynomial[DoubleField, BitPackedExponents]
		for i := 0; i < numPoly; i++ {
			var terms []Term[DoubleField, BitPackedExponents]
			for j := 0; j < numTerms; j++ {
				coeff := DoubleField{Value: rand.NextDouble()}
				exps := NewBitPackedExponents([6]uint8{uint8(rand.NextInt() % 4), uint8(rand.NextInt() % 4), uint8(rand.NextInt() % 4), 0, 0, 0})
				terms = append(terms, Term[DoubleField, BitPackedExponents]{Coefficient: coeff, Exponents: exps})
			}
			polys = append(polys, NewPolynomial(terms, order))
		}
		basis := NaiveGrobnerBasis(polys, order)
		fmt.Println("DoubleField + BitPackedExponents:")
		for i, poly := range basis {
			fmt.Printf("G%d: %v\n", i, poly)
		}
	case coeffType == 1 && expType == 0:
		// IntModP + VecExponents
		var polys []Polynomial[IntModP, VecExponents]
		for i := 0; i < numPoly; i++ {
			var terms []Term[IntModP, VecExponents]
			for j := 0; j < numTerms; j++ {
				coeff := NewIntModP(uint64(rand.NextInt()%int(modulus)), uint64(modulus))
				exps := NewVecExponents([]uint32{uint32(rand.NextInt() % 4), uint32(rand.NextInt() % 4), uint32(rand.NextInt() % 4)})
				terms = append(terms, Term[IntModP, VecExponents]{Coefficient: coeff, Exponents: exps})
			}
			polys = append(polys, NewPolynomial(terms, order))
		}
		basis := NaiveGrobnerBasis(polys, order)
		fmt.Println("IntModP + VecExponents:")
		for i, poly := range basis {
			fmt.Printf("G%d: %v\n", i, poly)
		}
	case coeffType == 1 && expType == 1:
		// IntModP + BitPackedExponents
		var polys []Polynomial[IntModP, BitPackedExponents]
		for i := 0; i < numPoly; i++ {
			var terms []Term[IntModP, BitPackedExponents]
			for j := 0; j < numTerms; j++ {
				coeff := NewIntModP(uint64(rand.NextInt()%int(modulus)), uint64(modulus))
				exps := NewBitPackedExponents([6]uint8{uint8(rand.NextInt() % 4), uint8(rand.NextInt() % 4), uint8(rand.NextInt() % 4), 0, 0, 0})
				terms = append(terms, Term[IntModP, BitPackedExponents]{Coefficient: coeff, Exponents: exps})
			}
			polys = append(polys, NewPolynomial(terms, order))
		}
		basis := NaiveGrobnerBasis(polys, order)
		fmt.Println("IntModP + BitPackedExponents:")
		for i, poly := range basis {
			fmt.Printf("G%d: %v\n", i, poly)
		}
	case coeffType == 2 && expType == 0:
		// SingleField + VecExponents
		var polys []Polynomial[SingleField, VecExponents]
		for i := 0; i < numPoly; i++ {
			var terms []Term[SingleField, VecExponents]
			for j := 0; j < numTerms; j++ {
				coeff := SingleField{Value: float32(rand.NextDouble())}
				exps := NewVecExponents([]uint32{uint32(rand.NextInt() % 4), uint32(rand.NextInt() % 4), uint32(rand.NextInt() % 4)})
				terms = append(terms, Term[SingleField, VecExponents]{Coefficient: coeff, Exponents: exps})
			}
			polys = append(polys, NewPolynomial(terms, order))
		}
		basis := NaiveGrobnerBasis(polys, order)
		fmt.Println("SingleField + VecExponents:")
		for i, poly := range basis {
			fmt.Printf("G%d: %v\n", i, poly)
		}
	case coeffType == 2 && expType == 1:
		// SingleField + BitPackedExponents
		var polys []Polynomial[SingleField, BitPackedExponents]
		for i := 0; i < numPoly; i++ {
			var terms []Term[SingleField, BitPackedExponents]
			for j := 0; j < numTerms; j++ {
				coeff := SingleField{Value: float32(rand.NextDouble())}
				exps := NewBitPackedExponents([6]uint8{uint8(rand.NextInt() % 4), uint8(rand.NextInt() % 4), uint8(rand.NextInt() % 4), 0, 0, 0})
				terms = append(terms, Term[SingleField, BitPackedExponents]{Coefficient: coeff, Exponents: exps})
			}
			polys = append(polys, NewPolynomial(terms, order))
		}
		basis := NaiveGrobnerBasis(polys, order)
		fmt.Println("SingleField + BitPackedExponents:")
		for i, poly := range basis {
			fmt.Printf("G%d: %v\n", i, poly)
		}
	default:
		fmt.Println("Invalid coefficient or exponent type selection.")
	}
}
