package specialized;

import java.util.*;

import helpers.LCG;
import helpers.ModInverse;


public class FiniteGrobner {
    public static TermOrder termOrder = TermOrder.Lex;
    public static int modulus;
    public static class Term {
        public int coefficient;
        public List<Integer> exponents; // Exponents for each variable

        public Term(int coefficient, List<Integer> exponents) {
            this.coefficient = ((coefficient % modulus) + modulus) % modulus;
            this.exponents = new ArrayList<>(exponents);
        }

        public String toString() {
            return String.format("%d * %s", coefficient, exponents.toString());
        }

        public int degree() {
            return exponents.stream().mapToInt(Integer::intValue).sum();
        }

        public int compareTo(Term other) {
            TermOrder order = termOrder;
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

        public Polynomial(List<Term> terms) {
            this.terms = new ArrayList<>(terms);
            this.terms.removeIf(t -> t.coefficient == 0);
            this.terms.sort((a, b) -> b.compareTo(a));
        }

        public Polynomial deepCopy() {
            List<Term> newTerms = new ArrayList<>();
            for (Term t : this.terms) {
                newTerms.add(new Term(t.coefficient, t.exponents));
            }
            return new Polynomial(newTerms);
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

        public Polynomial add(Polynomial other) {
            List<Term> result = this.deepCopy().terms;
            for (Term term : other.terms) {
                boolean found = false;
                for (Term resTerm : result) {
                    if (resTerm.exponents.equals(term.exponents)) {
                        resTerm.coefficient = (resTerm.coefficient + term.coefficient) % modulus;
                        found = true;
                        break;
                    }
                }
                if (!found) {
                    result.add(new Term(term.coefficient, term.exponents));
                }
            }
            return new Polynomial(result);
        }

        public Polynomial subtract(Polynomial other) {
            List<Term> result = this.deepCopy().terms;
            for (Term term : other.terms) {
                boolean found = false;
                for (Term resTerm : result) {
                    if (resTerm.exponents.equals(term.exponents)) {
                        resTerm.coefficient = (modulus + resTerm.coefficient - term.coefficient) % modulus;
                        found = true;
                        break;
                    }
                }
                if (!found) {
                    int negCoeff = (modulus - term.coefficient) % modulus;
                    result.add(new Term(negCoeff, term.exponents));
                }
            }
            return new Polynomial(result);
        }

        public Polynomial makeMonic() {
            if (terms.isEmpty()) return this.deepCopy();
            int leadCoeff = terms.get(0).coefficient;
            int inv = ModInverse.modInverse(leadCoeff, modulus);
            List<Term> newTerms = new ArrayList<>();
            for (Term t : terms) {
                int newCoeff = (t.coefficient * inv) % modulus;
                if (newCoeff < 0) newCoeff += modulus;
                newTerms.add(new Term(newCoeff, t.exponents));
            }
            return new Polynomial(newTerms);
        }

        public Polynomial multiplyByTerm(Term term) {
            List<Term> result = new ArrayList<>();
            for (Term t : terms) {
                List<Integer> newExponents = new ArrayList<>();
                for (int i = 0; i < t.exponents.size(); i++) {
                    newExponents.add(t.exponents.get(i) + term.exponents.get(i));
                }
                result.add(new Term((t.coefficient * term.coefficient) % modulus, newExponents));
            }
            return new Polynomial(result);
        }

        public Polynomial reduce(List<Polynomial> divisors) {
            // deep copy of self
            //Polynomial result = this.deepCopy(order);
            Polynomial result = new Polynomial(this.terms);
            List<Term> remainder = new ArrayList<>();
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
                        int coefficient = (leadingTerm.coefficient * ModInverse.modInverse(divisorLeadingTerm.coefficient, modulus)) % modulus;
                        List<Integer> exps = new ArrayList<>();
                        for (int i = 0; i < leadingTerm.exponents.size(); i++) {
                            exps.add(leadingTerm.exponents.get(i) - divisorLeadingTerm.exponents.get(i));
                        }
                        Term reductionTerm = new Term(coefficient, exps);
                        Polynomial scaledDivisor = divisor.multiplyByTerm(reductionTerm);
                        result = result.subtract(scaledDivisor);
                        reduced = true;
                        break;
                    }
                }
                if (!reduced) {
                    if (result.terms.isEmpty()) break;
                    //System.out.println("No reduction possible for leading term: " + result.terms.get(0).toString() + ", moving to remainder");
                    remainder.add(result.terms.get(0));
                    result.terms.remove(0);
                }
            }
            result.terms.addAll(remainder);
            return new Polynomial(result.terms);
        }

        public static Polynomial sPolynomial(Polynomial p1, Polynomial p2) {
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
                scaledTermsP1.add(new Term(t.coefficient, newExponents));
            }
            List<Term> scaledTermsP2 = new ArrayList<>();
            for (Term t : p2.terms) {
                List<Integer> newExponents = new ArrayList<>();
                for (int i = 0; i < n; i++) {
                    newExponents.add(t.exponents.get(i) + scaleFactorP2.get(i));
                }
                scaledTermsP2.add(new Term(t.coefficient, newExponents));
            }
            Polynomial scaledP1 = new Polynomial(scaledTermsP1);
            Polynomial scaledP2 = new Polynomial(scaledTermsP2);
            return scaledP1.subtract(scaledP2);
        }
    }

    public static List<Polynomial> naiveGrobnerBasis(List<Polynomial> polynomials) {
        List<Polynomial> basis = new ArrayList<>();
        Set<Polynomial> basisSet = new HashSet<>();
        for (Polynomial poly : polynomials) {
            //basis.add(poly.deepCopy(order));
            basis.add(poly);
        }
        List<int[]> pairs = new ArrayList<>();
        for (int i = 0; i < basis.size(); i++) {
            for (int j = i + 1; j < basis.size(); j++) {
                pairs.add(new int[]{i, j});
            }
        }
        while (!pairs.isEmpty()) {
            int[] pair = pairs.remove(0);
            int i = pair[0]; int j = pair[1];
            Polynomial sPoly = Polynomial.sPolynomial(basis.get(i), basis.get(j));
            Polynomial reduced = sPoly.reduce(basis);
            // If non-trivial and new, add to basis and enqueue new pairs
            if (!reduced.terms.isEmpty() && !basisSet.contains(reduced)) {
                basisSet.add(reduced);
                int newPolyIdx = basis.size();
                basis.add(reduced);
                
                // Add pairs between new polynomial and all existing ones
                for (int k = 0; k < newPolyIdx; k++) {
                    pairs.add(new int[]{k, newPolyIdx});
                }
            }
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
            Polynomial reduced = poly.reduce(basisExcludingSelf);
            if (!reduced.terms.isEmpty() && !reducedBasis.contains(reduced)) {
                reducedBasis.add(reduced.makeMonic());
            }
        }
        return reducedBasis;
    }

    public static boolean areBasesEquivalent(List<Polynomial> setA, List<Polynomial> setB, TermOrder order) {
        for (Polynomial poly : setA) {
            Polynomial reduced = poly.reduce(setB);
            if (!reduced.terms.isEmpty()) return false;
        }
        for (Polynomial poly : setB) {
            Polynomial reduced = poly.reduce(setA);
            if (!reduced.terms.isEmpty()) return false;
        }
        return true;
    }

    public static void main(String[] args) {
        // let mode = 0 be for testing
        int mode = 0;
        FiniteGrobner.termOrder = TermOrder.Lex;
        if (args.length > 1) {
            int orderArg = Integer.parseInt(args[1]);
            switch (orderArg) {
                case 0: FiniteGrobner.termOrder = TermOrder.Lex; break;
                case 1: FiniteGrobner.termOrder = TermOrder.GrLex; break;
                case 2: FiniteGrobner.termOrder = TermOrder.RevLex; break;
                default: FiniteGrobner.termOrder = TermOrder.Lex;
            }
        }
        if (mode != 0) {
            int numPolynomials = args.length > 0 ? Integer.parseInt(args[0]) : 3;
            int modulus = 13;
            FiniteGrobner.modulus = modulus;
            LCG rand = new LCG(12345, 1345, 16645, 1013904);
            List<Polynomial> inputBasis = new ArrayList<>();
            int numTerms = 3;
            for (int i = 0; i < numPolynomials; i++) {
                List<Term> terms = new ArrayList<>();
                for (int j = 0; j < numTerms; j++) {
                    int coefficient = rand.nextInt() % 9 + 1;
                    List<Integer> exponents = Arrays.asList(rand.nextInt() % 4, rand.nextInt() % 4, rand.nextInt() % 4);
                    terms.add(new Term(coefficient, exponents));
                }
                inputBasis.add(new Polynomial(terms));
            }
            List<Polynomial> basis = naiveGrobnerBasis(inputBasis);
            System.out.println(basis.size());
            
            System.out.println("Computed Grobner Basis Polynomials:");
            for (Polynomial poly : basis) {
                System.out.println(poly.terms);
            }
            
        } else {
            int modulus = 7;
            FiniteGrobner.modulus = modulus;
            
            int n = args.length > 0 ? Integer.parseInt(args[0]) : 5;
            if (n == 4) {
                System.out.println("Java Specialized finite coeff vector exp cyclic 4");
                // x + y + z + w
                List<Term> terms1 = new ArrayList<>();
                terms1.add(new Term(1, Arrays.asList(1, 0, 0, 0)));
                terms1.add(new Term(1, Arrays.asList(0, 1, 0, 0)));
                terms1.add(new Term(1, Arrays.asList(0, 0, 1, 0)));
                terms1.add(new Term(1, Arrays.asList(0, 0, 0, 1)));
                Polynomial p1 = new Polynomial(terms1);

                // x*y + y*z + z*w + w*x
                List<Term> terms2 = new ArrayList<>();
                terms2.add(new Term(1, Arrays.asList(1, 1, 0, 0)));
                terms2.add(new Term(1, Arrays.asList(0, 1, 1, 0)));
                terms2.add(new Term(1, Arrays.asList(0, 0, 1, 1)));
                terms2.add(new Term(1, Arrays.asList(1, 0, 0, 1)));
                Polynomial p2 = new Polynomial(terms2);

                // x*y*z + y*z*w + z*w*x + w*x*y
                List<Term> terms3 = new ArrayList<>();
                terms3.add(new Term(1, Arrays.asList(1, 1, 1, 0)));
                terms3.add(new Term(1, Arrays.asList(0, 1, 1, 1)));
                terms3.add(new Term(1, Arrays.asList(1, 0, 1, 1)));
                terms3.add(new Term(1, Arrays.asList(1, 1, 0, 1)));
                Polynomial p3 = new Polynomial(terms3);

                // x*y*z*w - 1
                List<Term> terms4 = new ArrayList<>();
                terms4.add(new Term(1, Arrays.asList(1, 1, 1, 1)));
                terms4.add(new Term(modulus-1, Arrays.asList(0, 0, 0, 0)));
                Polynomial p4 = new Polynomial(terms4);


                for (int i = 0; i < 10; i++) {
                    List<Polynomial> inputBasis = Arrays.asList(p1, p2, p3, p4);
                    List<Polynomial> basis = naiveGrobnerBasis(inputBasis);
                    System.out.println("Iteration " + i + " complete");
                    if (i == 9) {
                        System.out.println("Cyclic 4 Grobner Basis:");
                        for (Polynomial poly : basis) {
                            System.out.println(poly.terms);
                        }
                    }
                }
            }
            else if (n == 5) {
                System.out.println("Java Specialized finite coeff vector exp cyclic 5");
                // x + y + z + w + v
                List<Term> terms1 = new ArrayList<>();
                terms1.add(new Term(1, Arrays.asList(1, 0, 0, 0, 0)));
                terms1.add(new Term(1, Arrays.asList(0, 1, 0, 0, 0)));
                terms1.add(new Term(1, Arrays.asList(0, 0, 1, 0, 0)));
                terms1.add(new Term(1, Arrays.asList(0, 0, 0, 1, 0)));
                terms1.add(new Term(1, Arrays.asList(0, 0, 0, 0, 1)));
                Polynomial p1 = new Polynomial(terms1);

                // xy + yz + zw + wv + vx
                List<Term> terms2 = new ArrayList<>();
                terms2.add(new Term(1, Arrays.asList(1, 1, 0, 0, 0)));
                terms2.add(new Term(1, Arrays.asList(0, 1, 1, 0, 0)));
                terms2.add(new Term(1, Arrays.asList(0, 0, 1, 1, 0)));
                terms2.add(new Term(1, Arrays.asList(0, 0, 0, 1, 1)));
                terms2.add(new Term(1, Arrays.asList(1, 0, 0, 0, 1)));
                Polynomial p2 = new Polynomial(terms2);

                // xyz + yzw + zwv + wvx + vxy
                List<Term> terms3 = new ArrayList<>();
                terms3.add(new Term(1, Arrays.asList(1, 1, 1, 0, 0)));
                terms3.add(new Term(1, Arrays.asList(0, 1, 1, 1, 0)));
                terms3.add(new Term(1, Arrays.asList(0, 0, 1, 1, 1)));
                terms3.add(new Term(1, Arrays.asList(1, 0, 0, 1, 1)));
                terms3.add(new Term(1, Arrays.asList(1, 1, 0, 0, 1)));
                Polynomial p3 = new Polynomial(terms3);

                // xyzw + yzwv + zwvx + wvxy + vxyz
                List<Term> terms4 = new ArrayList<>();
                terms4.add(new Term(1, Arrays.asList(1, 1, 1, 1, 0)));
                terms4.add(new Term(1, Arrays.asList(0, 1, 1, 1, 1)));
                terms4.add(new Term(1, Arrays.asList(1, 0, 1, 1, 1)));
                terms4.add(new Term(1, Arrays.asList(1, 1, 0, 1, 1)));
                terms4.add(new Term(1, Arrays.asList(1, 1, 1, 0, 1)));
                Polynomial p4 = new Polynomial(terms4);

                // xyzwv - 1
                List<Term> terms5 = new ArrayList<>();
                terms5.add(new Term(1, Arrays.asList(1, 1, 1, 1, 1)));
                terms5.add(new Term(modulus-1, Arrays.asList(0, 0, 0, 0, 0)));
                Polynomial p5 = new Polynomial(terms5);

                for (int i = 0; i < 10; i++) {
                    List<Polynomial> inputBasis = Arrays.asList(p1, p2, p3, p4, p5);
                    List<Polynomial> basis = naiveGrobnerBasis(inputBasis);
                    System.out.println("Iteration " + i + " complete");
                    if (i == 9) {
                        System.out.println("Cyclic 5 Grobner Basis:");
                        for (Polynomial poly : basis) {
                            System.out.println(poly.terms);
                        }
                    }
                }
            }
            else if (n == 6) {
                System.out.println("Java Specialized finite coeff vector exp cyclic 6");
                // x + y + z + w + v + u
                List<Term> terms1 = new ArrayList<>();
                terms1.add(new Term(1, Arrays.asList(1, 0, 0, 0, 0, 0)));
                terms1.add(new Term(1, Arrays.asList(0, 1, 0, 0, 0, 0)));
                terms1.add(new Term(1, Arrays.asList(0, 0, 1, 0, 0, 0)));
                terms1.add(new Term(1, Arrays.asList(0, 0, 0, 1, 0, 0)));
                terms1.add(new Term(1, Arrays.asList(0, 0, 0, 0, 1, 0)));
                terms1.add(new Term(1, Arrays.asList(0, 0, 0, 0, 0, 1)));
                Polynomial p1 = new Polynomial(terms1);

                // xy + yz + zw + wv + vu + ux
                List<Term> terms2 = new ArrayList<>();
                terms2.add(new Term(1, Arrays.asList(1, 1, 0, 0, 0, 0)));
                terms2.add(new Term(1, Arrays.asList(0, 1, 1, 0, 0, 0)));
                terms2.add(new Term(1, Arrays.asList(0, 0, 1, 1, 0, 0)));
                terms2.add(new Term(1, Arrays.asList(0, 0, 0, 1, 1, 0)));
                terms2.add(new Term(1, Arrays.asList(0, 0, 0, 0, 1, 1)));
                terms2.add(new Term(1, Arrays.asList(1, 0, 0, 0, 0, 1)));
                Polynomial p2 = new Polynomial(terms2);

                // xyz + yzw + zwv + wvu + vux + uxy
                List<Term> terms3 = new ArrayList<>();
                terms3.add(new Term(1, Arrays.asList(1, 1, 1, 0, 0, 0)));
                terms3.add(new Term(1, Arrays.asList(0, 1, 1, 1, 0, 0)));
                terms3.add(new Term(1, Arrays.asList(0, 0, 1, 1, 1, 0)));
                terms3.add(new Term(1, Arrays.asList(0, 0, 0, 1, 1, 1)));
                terms3.add(new Term(1, Arrays.asList(1, 0, 0, 0, 1, 1)));
                terms3.add(new Term(1, Arrays.asList(1, 1, 0, 0, 0, 1)));
                Polynomial p3 = new Polynomial(terms3);

                // xyzw + yzwv + zwvu + wvux + vuxy + uxyz
                List<Term> terms4 = new ArrayList<>();
                terms4.add(new Term(1, Arrays.asList(1, 1, 1, 1, 0, 0)));
                terms4.add(new Term(1, Arrays.asList(0, 1, 1, 1, 1, 0)));
                terms4.add(new Term(1, Arrays.asList(0, 0, 1, 1, 1, 1)));
                terms4.add(new Term(1, Arrays.asList(1, 0, 0, 1, 1, 1)));
                terms4.add(new Term(1, Arrays.asList(1, 1, 0, 0, 1, 1)));
                terms4.add(new Term(1, Arrays.asList(1, 1, 1, 0, 0, 1)));
                Polynomial p4 = new Polynomial(terms4);

                // xyzwv + yzwvu + zwvux + wvuxy + vuxyz + uxyxzw
                List<Term> terms5 = new ArrayList<>();
                terms5.add(new Term(1, Arrays.asList(1, 1, 1, 1, 1, 0)));
                terms5.add(new Term(1, Arrays.asList(0, 1, 1, 1, 1, 1)));
                terms5.add(new Term(1, Arrays.asList(1, 0, 1, 1, 1, 1)));
                terms5.add(new Term(1, Arrays.asList(1, 1, 0, 1, 1, 1)));
                terms5.add(new Term(1, Arrays.asList(1, 1, 1, 0, 1, 1)));
                terms5.add(new Term(1, Arrays.asList(1, 1, 1, 1, 0, 1)));
                Polynomial p5 = new Polynomial(terms5);

                // xyzwvu - 1
                List<Term> terms6 = new ArrayList<>();
                terms6.add(new Term(1, Arrays.asList(1, 1, 1, 1, 1, 1)));
                terms6.add(new Term(modulus-1, Arrays.asList(0, 0, 0, 0, 0, 0)));
                Polynomial p6 = new Polynomial(terms6);

                for (int i = 0; i < 10; i++) {
                    List<Polynomial> inputBasis = Arrays.asList(p1, p2, p3, p4, p5, p6);
                    List<Polynomial> basis = naiveGrobnerBasis(inputBasis);
                    System.out.println("Iteration " + i + " complete");
                    if (i == 9) {
                        System.out.println("Cyclic 6 Grobner Basis:");
                        for (Polynomial poly : basis) {
                            System.out.println(poly.terms);
                        }
                    }
                }
            }
            
        }
    }
}
