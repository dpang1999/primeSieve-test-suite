package specialized;

import java.util.*;

import helpers.LCG;
import helpers.ModInverse;


public class FiniteGrobner {
    public static class Term {
        public int coefficient;
        public List<Integer> exponents; // Exponents for each variable
        public int modulus;

        public Term(int coefficient, List<Integer> exponents, int modulus) {
            this.coefficient = coefficient;
            this.exponents = new ArrayList<>(exponents);
            this.modulus = modulus;
        }

        public String toString() {
            return String.format("%d * %s", coefficient, exponents.toString());
        }

        public int degree() {
            return exponents.stream().mapToInt(Integer::intValue).sum();
        }

        public int compareTo(Term other, TermOrder order) {
            switch (order) {
                case Lex:
                    for (int i = 0; i < exponents.size(); i++) {
                        int cmp = Integer.compare(exponents.get(i), other.exponents.get(i));
                        if (cmp != 0) return cmp;
                    }
                    return 0;
                case GrLex:
                    int selfDegree = this.degree();
                    int otherDegree = other.degree();
                    int cmp = Integer.compare(selfDegree, otherDegree);
                    if (cmp != 0) return cmp;
                    for (int i = 0; i < exponents.size(); i++) {
                        cmp = Integer.compare(exponents.get(i), other.exponents.get(i));
                        if (cmp != 0) return cmp;
                    }
                    return 0;
                case RevLex:
                    selfDegree = this.degree();
                    otherDegree = other.degree();
                    cmp = Integer.compare(selfDegree, otherDegree);
                    if (cmp != 0) return cmp;
                    for (int i = exponents.size() - 1; i >= 0; i--) {
                        cmp = Integer.compare(other.exponents.get(i), exponents.get(i));
                        if (cmp != 0) return cmp;
                    }
                    return 0;
                default:
                    return 0;
            }
        }

        @Override
        public boolean equals(Object o) {
            if (this == o) return true;
            if (!(o instanceof Term)) return false;
            Term term = (Term) o;
            return coefficient == term.coefficient && Objects.equals(exponents, term.exponents);
        }

        @Override
        public int hashCode() {
            return Objects.hash(coefficient, exponents);
        }
    }

    public enum TermOrder {
        Lex, GrLex, RevLex
    }

    public static class Polynomial {
        public List<Term> terms;

        public Polynomial(List<Term> terms, TermOrder order) {
            this.terms = new ArrayList<>(terms);
            this.terms.removeIf(t -> t.coefficient == 0);
            this.terms.sort((a, b) -> b.compareTo(a, order));
        }

        public Polynomial deepCopy(TermOrder order) {
            List<Term> newTerms = new ArrayList<>();
            for (Term t : this.terms) {
                newTerms.add(new Term(t.coefficient, t.exponents, t.modulus));
            }
            return new Polynomial(newTerms, order);
        }

        @Override
        public String toString() {
            if (terms.isEmpty()) return "0";
            StringBuilder sb = new StringBuilder();
            for (int i = 0; i < terms.size(); i++) {
                if (i > 0) sb.append(" + ");
                sb.append(terms.get(i).toString());
            }
            return sb.toString();
        }

        public Polynomial add(Polynomial other, TermOrder order) {
            List<Term> result = this.deepCopy(order).terms;
            for (Term term : other.terms) {
                boolean found = false;
                for (Term resTerm : result) {
                    if (resTerm.exponents.equals(term.exponents)) {
                        resTerm.coefficient = (resTerm.coefficient + term.coefficient) % term.modulus;
                        found = true;
                        break;
                    }
                }
                if (!found) {
                    result.add(new Term(term.coefficient, term.exponents, term.modulus));
                }
            }
            return new Polynomial(result, order);
        }

        public Polynomial subtract(Polynomial other, TermOrder order) {
            List<Term> result = this.deepCopy(order).terms;
            for (Term term : other.terms) {
                boolean found = false;
                for (Term resTerm : result) {
                    if (resTerm.exponents.equals(term.exponents)) {
                        resTerm.coefficient = (resTerm.coefficient - term.coefficient) % term.modulus;
                        found = true;
                        break;
                    }
                }
                if (!found) {
                    int negCoeff = (term.modulus - (term.coefficient % term.modulus)) % term.modulus;
                    result.add(new Term(negCoeff, term.exponents, term.modulus));
                }
            }
            return new Polynomial(result, order);
        }

        public Polynomial multiplyByTerm(Term term, TermOrder order) {
            List<Term> result = new ArrayList<>();
            for (Term t : terms) {
                List<Integer> newExponents = new ArrayList<>();
                for (int i = 0; i < t.exponents.size(); i++) {
                    newExponents.add(t.exponents.get(i) + term.exponents.get(i));
                }
                result.add(new Term((t.coefficient * term.coefficient) % term.modulus, newExponents, term.modulus));
            }
            return new Polynomial(result, order);
        }

        public Polynomial reduce(List<Polynomial> divisors, TermOrder order) {
            // deep copy of self
            //Polynomial result = this.deepCopy(order);
            Polynomial result = new Polynomial(this.terms, order);
            while (true) {
                boolean reduced = false;
                for (Polynomial divisor : divisors) {
                    if (result.terms.isEmpty() || divisor.terms.isEmpty()) continue;
                    Term leadingTerm = result.terms.get(0);
                    Term divisorLeadingTerm = divisor.terms.get(0);
                    boolean canReduce = true;
                    for (int i = 0; i < leadingTerm.exponents.size(); i++) {
                        if (leadingTerm.exponents.get(i) < divisorLeadingTerm.exponents.get(i)) {
                            canReduce = false;
                            break;
                        }
                    }
                    if (canReduce) {
                        int modulus = leadingTerm.modulus;
                        int coefficient = (leadingTerm.coefficient * ModInverse.modInverse(divisorLeadingTerm.coefficient, modulus)) % modulus;
                        List<Integer> exps = new ArrayList<>();
                        for (int i = 0; i < leadingTerm.exponents.size(); i++) {
                            exps.add(leadingTerm.exponents.get(i) - divisorLeadingTerm.exponents.get(i));
                        }
                        Term reductionTerm = new Term(coefficient, exps, leadingTerm.modulus);
                        Polynomial scaledDivisor = divisor.multiplyByTerm(reductionTerm, order);
                        result = result.subtract(scaledDivisor, order);
                        reduced = true;
                        break;
                    }
                }
                if (!reduced) break;
            }
            return new Polynomial(result.terms, order);
        }

        public static Polynomial sPolynomial(Polynomial p1, Polynomial p2, TermOrder order) {
            Term leadingTermP1 = p1.terms.get(0);
            Term leadingTermP2 = p2.terms.get(0);
            int n = leadingTermP1.exponents.size();
            List<Integer> lcmExponents = new ArrayList<>();
            for (int i = 0; i < n; i++) {
                lcmExponents.add(Math.max(leadingTermP1.exponents.get(i), leadingTermP2.exponents.get(i)));
            }
            List<Integer> scaleFactorP1 = new ArrayList<>();
            List<Integer> scaleFactorP2 = new ArrayList<>();
            for (int i = 0; i < n; i++) {
                scaleFactorP1.add(lcmExponents.get(i) - leadingTermP1.exponents.get(i));
                scaleFactorP2.add(lcmExponents.get(i) - leadingTermP2.exponents.get(i));
            }
            // Scale all terms in p1 and p2 by the scale factor
            List<Term> scaledTermsP1 = new ArrayList<>();
            for (Term t : p1.terms) {
                List<Integer> newExponents = new ArrayList<>();
                for (int i = 0; i < n; i++) {
                    newExponents.add(t.exponents.get(i) + scaleFactorP1.get(i));
                }
                scaledTermsP1.add(new Term(t.coefficient, newExponents, t.modulus));
            }
            List<Term> scaledTermsP2 = new ArrayList<>();
            for (Term t : p2.terms) {
                List<Integer> newExponents = new ArrayList<>();
                for (int i = 0; i < n; i++) {
                    newExponents.add(t.exponents.get(i) + scaleFactorP2.get(i));
                }
                scaledTermsP2.add(new Term(t.coefficient, newExponents, t.modulus));
            }
            Polynomial scaledP1 = new Polynomial(scaledTermsP1, order);
            Polynomial scaledP2 = new Polynomial(scaledTermsP2, order);
            return scaledP1.subtract(scaledP2, order);
        }
    }

    public static List<Polynomial> naiveGrobnerBasis(List<Polynomial> polynomials, TermOrder order) {
        List<Polynomial> basis = new ArrayList<>();
        for (Polynomial poly : polynomials) {
            //basis.add(poly.deepCopy(order));
            basis.add(poly);
        }
        Set<Polynomial> basisSet = new HashSet<>(basis);
        while (true) {
            boolean added = false;
            int basisLen = basis.size();
            for (int i = 0; i < basisLen; i++) {
                for (int j = i + 1; j < basisLen; j++) {
                    Polynomial sPoly = Polynomial.sPolynomial(basis.get(i), basis.get(j), order);
                    // print basis terms and s polynomial
                    //System.out.print("Basis 1: " + basis.get(i).terms);
                    //System.out.print(" | Basis 2: " + basis.get(j).terms);
                    //System.out.println(" | S-Polynomial: " + sPoly.terms);
                    Polynomial reduced = sPoly.reduce(basis, order);
                    if (!reduced.terms.isEmpty() && !basisSet.contains(reduced)) {
                        //basis.add(reduced.deepCopy(order));
                        //basisSet.add(reduced.deepCopy(order));
                        basis.add(reduced);
                        basisSet.add(reduced);
                        added = true;
                    }
                }
            }
            if (!added) break;
        }
        List<Polynomial> reducedBasis = new ArrayList<>();
         for (Polynomial poly : basis) {
            //System.out.println("Basis: " + basis.toString());
            List<Polynomial> basisExcludingSelf = new ArrayList<>();
            for (Polynomial p : basis) {
                if (!p.equals(poly)) {
                    //basisExcludingSelf.add(p.deepCopy(order));
                    basisExcludingSelf.add(p);
                }
            }
            Polynomial reduced = poly.reduce(basisExcludingSelf, order);
            if (!reduced.terms.isEmpty() && !reducedBasis.contains(reduced)) {
                reducedBasis.add(reduced);
            }
        }
        return reducedBasis;
    }

    public static boolean areBasesEquivalent(List<Polynomial> setA, List<Polynomial> setB, TermOrder order) {
        for (Polynomial poly : setA) {
            Polynomial reduced = poly.reduce(setB, order);
            if (!reduced.terms.isEmpty()) return false;
        }
        for (Polynomial poly : setB) {
            Polynomial reduced = poly.reduce(setA, order);
            if (!reduced.terms.isEmpty()) return false;
        }
        return true;
    }

    public static void main(String[] args) {
        // let mode = 0 be for testing
        int mode = 1;
        TermOrder order = TermOrder.Lex;
        if (args.length > 1) {
            int orderArg = Integer.parseInt(args[1]);
            switch (orderArg) {
                case 0: order = TermOrder.Lex; break;
                case 1: order = TermOrder.GrLex; break;
                case 2: order = TermOrder.RevLex; break;
                default: order = TermOrder.Lex;
            }
        }
        if (mode != 0) {
            int numPolynomials = args.length > 0 ? Integer.parseInt(args[0]) : 3;
            int modulus = 13;
            LCG rand = new LCG(12345, 1345, 16645, 1013904);
            List<Polynomial> inputBasis = new ArrayList<>();
            int numTerms = 3;
            for (int i = 0; i < numPolynomials; i++) {
                List<Term> terms = new ArrayList<>();
                for (int j = 0; j < numTerms; j++) {
                    int coefficient = rand.nextInt() % 9 + 1;
                    List<Integer> exponents = Arrays.asList(rand.nextInt() % 4, rand.nextInt() % 4, rand.nextInt() % 4);
                    terms.add(new Term(coefficient, exponents, modulus));
                }
                inputBasis.add(new Polynomial(terms, order));
            }
            List<Polynomial> basis = naiveGrobnerBasis(inputBasis, order);
            System.out.println(basis.size());
            
            System.out.println("Computed Grobner Basis Polynomials:");
            for (Polynomial poly : basis) {
                System.out.println(poly.terms);
            }
            
        } else {
            int modulus = 13;
            // x^2 - y
            List<Term> terms1 = new ArrayList<>();
            terms1.add(new Term(1, Arrays.asList(2, 0), modulus));
            terms1.add(new Term(-1, Arrays.asList(0, 1), modulus));
            Polynomial p1 = new Polynomial(terms1, order);
            // xy - 1   
            List<Term> terms2 = new ArrayList<>();
            terms2.add(new Term(1, Arrays.asList(1, 1), modulus));
            terms2.add(new Term(-1, Arrays.asList(0, 0), modulus));
            Polynomial p2 = new Polynomial(terms2, order);
            List<Polynomial> inputBasis = Arrays.asList(p1, p2);
            List<Polynomial> basis = naiveGrobnerBasis(inputBasis, order);
            System.out.println("Computed Grobner Basis Polynomials:");
            for (Polynomial poly : basis) {
                System.out.println(poly.terms);
            }


            // x^3 + y^3 + z^3
            List<Term> terms3 = new ArrayList<>();
            terms3.add(new Term(1, Arrays.asList(3, 0, 0), modulus));
            terms3.add(new Term(1, Arrays.asList(0, 3, 0), modulus));
            terms3.add(new Term(1, Arrays.asList(0, 0, 3), modulus));
            Polynomial p3 = new Polynomial(terms3, order);
            // x*y + y*z + z*x
            List<Term> terms4 = new ArrayList<>();
            terms4.add(new Term(1, Arrays.asList(1, 1, 0), modulus));
            terms4.add(new Term(1, Arrays.asList(0, 1, 1), modulus));
            terms4.add(new Term(1, Arrays.asList(1, 0, 1), modulus));
            Polynomial p4 = new Polynomial(terms4, order);
            // x+y+z
            List<Term> terms5 = new ArrayList<>();
            terms5.add(new Term(1, Arrays.asList(1, 0, 0), modulus));
            terms5.add(new Term(1, Arrays.asList(0, 1, 0), modulus));
            terms5.add(new Term(1, Arrays.asList(0, 0, 1), modulus));
            Polynomial p5 = new Polynomial(terms5   , order);
            List<Polynomial> inputBasis2 = Arrays.asList(p3, p4, p5);
            List<Polynomial> basis2 = naiveGrobnerBasis(inputBasis2, order);
            System.out.println("Computed Grobner Basis Polynomials for second example:");
            for (Polynomial poly : basis2) {
                System.out.println(poly.terms);
            }
        }
    }
}
