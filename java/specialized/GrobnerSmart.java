package specialized;

import java.util.*;

import helpers.LCG;
import helpers.ModInverse;

public class GrobnerSmart {
    public static TermOrder termOrder = TermOrder.Lex;
    public static int modulus;
    public static class Term {
        public long coefficient;
        public long exponents; // Bitpacked: [63..48]=degree, [47..40]=e0, [39..32]=e1, ... [7..0]=e5

        public Term(long coefficient, long exponents) {
            this.coefficient = coefficient;
            this.exponents = exponents;
        }

        // Create a term from array of exponents (length 6)
        public static Term fromExponents(long coefficient, int[] exps) {
            long packed = 0L;
            int degree = 0;
            for (int i = 0; i < 6; i++) {
                packed |= ((long) (exps[i] & 0xFF)) << (40 - 8 * i);
                degree += exps[i];
            }
            packed |= ((long) degree) << 48;
            return new Term(coefficient, packed);
        }

       
        // LCM of two exponent vectors
        public static long lcm(long a, long b) {
            long ea = a & 0x0000FFFFFFFFFFFFL;
            long eb = b & 0x0000FFFFFFFFFFFFL;
            long lcmExponents = 0;
            int degree = 0;
            for (int i = 0; i < 6; i++) {
                int shift = 40 - (8 * i);
                int expA = (int) ((ea >> shift) & 0xFF);
                int expB = (int) ((eb >> shift) & 0xFF);
                int maxExp = Math.max(expA, expB);
                lcmExponents |= ((long) maxExp) << shift;
                degree += maxExp;
            }
            return lcmExponents | (((long) degree) << 48);
        }

        // Compare terms by order
        public int compareTo(Term other) {
            TermOrder order = GrobnerSmart.termOrder;
            int cmp;
            int selfDeg = (int) ((this.exponents >> 48) & 0xFFFF);
            int otherDeg = (int) ((other.exponents >> 48) & 0xFFFF);
            switch (order) {
                case Lex:
                    //System.out.println("Comparing: " + hexExponents(this.exponents) + " to " + hexExponents(other.exponents));
                    cmp = Long.compare(this.exponents & 0x0000FFFFFFFFFFFFL, other.exponents & 0x0000FFFFFFFFFFFFL);
                    //System.out.println(cmp);
                    return cmp;
                case GrLex:
                    cmp = Integer.compare(selfDeg, otherDeg);
                    if (cmp != 0) return cmp;
                    // else degrees are equal 
                    cmp = Long.compare(this.exponents & 0x0000FFFFFFFFFFFFL, other.exponents & 0x0000FFFFFFFFFFFFL);
                    return cmp;
                case RevLex:
                    cmp = Integer.compare(selfDeg, otherDeg);
                    if (cmp != 0) return cmp;
                    // else degrees are equal
                    cmp = Long.compare(other.exponents & 0x0000FFFFFFFFFFFFL, this.exponents & 0x0000FFFFFFFFFFFFL);
                    return cmp;
                default:
                    return 0;
            }
        }

        // Can this term be reduced by divisor?
        public boolean canReduce(Term divisor) {
            long selfExps = this.exponents & 0x0000FFFFFFFFFFFFL;
            long divExps = divisor.exponents & 0x0000FFFFFFFFFFFFL;
            for (int i = 0; i < 6; i++) {
                if ((selfExps >> (i*8) & 0xFF) < (divExps >> (i*8) & 0xFF)) return false;
            }
            return true;
        }

        @Override
        public boolean equals(Object o) {
            if (this == o) return true;
            if (!(o instanceof Term)) return false;
            Term term = (Term) o;
            return coefficient == term.coefficient && exponents == term.exponents;
        }

        @Override
        public int hashCode() {
            return Objects.hash(coefficient, exponents);
        }

        @Override
        public String toString() {
            return String.format("%d * %s", coefficient, Arrays.toString(unpackExponents(exponents)));
        }

        // Print packed exponents in hex (ignoring degree)
        public static String hexExponents(long packed) {
            StringBuilder sb = new StringBuilder();
            for (int i = 0; i < 6; i++) {
                int shift = 40 - 8 * i;
                int exp = (int) ((packed >> shift) & 0xFF);
                sb.append(String.format("%02X ", exp));
            }
            return sb.toString().trim();
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
            //System.out.println("Before sort: " + this.terms);
            this.terms.sort((a, b) -> b.compareTo(a));
            //System.out.println("After sort: " + this.terms);
        }

        public Polynomial deepCopy() {
            List<Term> newTerms = new ArrayList<>();
            for (Term t : this.terms) {
                newTerms.add(new Term(t.coefficient, t.exponents));
            }
            return new Polynomial(newTerms);
        }

        @Override
        public boolean equals(Object o) {
            if (this == o) return true;
            if (!(o instanceof Polynomial)) return false;
            Polynomial poly = (Polynomial) o;
            return Objects.equals(terms, poly.terms);
        }

        @Override
        public int hashCode() {
            return Objects.hash(terms);
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
                    if (resTerm.exponents == term.exponents) {
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
                    if (resTerm.exponents == term.exponents) {
                        resTerm.coefficient = (modulus + resTerm.coefficient - term.coefficient) % modulus;
                        found = true;
                        break;
                    }
                }
                if (!found) {
                    result.add(new Term((modulus - term.coefficient) % modulus, term.exponents));
                }
            }
            return new Polynomial(result);
        }
        public Polynomial makeMonic() {
            if (terms.isEmpty()) return this.deepCopy();
            long leadCoeff = terms.get(0).coefficient;
            long inv = ModInverse.modInverse(leadCoeff, modulus);
            List<Term> newTerms = new ArrayList<>();
            for (Term t : terms) {
                long newCoeff = (t.coefficient * inv) % modulus;
                if (newCoeff < 0) newCoeff += modulus;
                newTerms.add(new Term(newCoeff, t.exponents));
            }
            return new Polynomial(newTerms);
        }

        public Polynomial multiplyByTerm(Term term) {
            List<Term> result = new ArrayList<>();
            for (Term t : this.terms) {
                long exp = t.exponents + term.exponents;
                long coefficient = (t.coefficient * term.coefficient) % modulus;
                result.add(new Term(coefficient, exp));
            }
            return new Polynomial(result);
        }

        public Polynomial reduce(List<Polynomial> divisors) {
            // deep copy of self
            Polynomial result = new Polynomial(this.terms);
            List<Term> remainder = new ArrayList<>();
            //Polynomial result = new Polynomial(this.terms, order);
            while (true) {
                boolean reduced = false;
                //System.out.println("Current polynomial to reduce: " + result.toString());
                for (Polynomial divisor : divisors) {
                    // debug
                    //System.out.println("Attempting to reduce by divisor: " + divisor.toString());
                    if (result.terms.isEmpty() || divisor.terms.isEmpty()) continue;
                    Term leadingTerm = result.terms.get(0);
                    Term divisorLeadingTerm = divisor.terms.get(0);
                    if (leadingTerm.canReduce(divisorLeadingTerm)) {
                        // debug
                        //System.out.println("leading term: " + leadingTerm.toString() + " can be reduced by divisor leading term: " + divisorLeadingTerm.toString());
                        long exponents = leadingTerm.exponents - divisorLeadingTerm.exponents;
                        long coefficient = (leadingTerm.coefficient * ModInverse.modInverse(divisorLeadingTerm.coefficient, (long)modulus)) % modulus;
                        Term reductionTerm = new Term(coefficient, exponents);
                        //System.out.println("reduction term: " + reductionTerm.toString());
                        Polynomial scaledDivisor = divisor.multiplyByTerm(reductionTerm);
                        //System.out.println("scaled divisor: " + scaledDivisor.toString());
                        result = result.subtract(scaledDivisor);
                        reduced = true;
                        //System.out.println("Reduced polynomial: " + result.toString());
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
            long lcmExponents = Term.lcm(leadingTermP1.exponents, leadingTermP2.exponents);
            long exp1 = lcmExponents - leadingTermP1.exponents;
            long exp2 = lcmExponents - leadingTermP2.exponents;
            Term multTerm1 = new Term(1L, exp1);
            Term multTerm2 = new Term(1L, exp2);
            Polynomial term1 = p1.multiplyByTerm(multTerm1);
            Polynomial term2 = p2.multiplyByTerm(multTerm2);
            return term1.subtract(term2);
            
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
        //System.out.println("Checkpoint");
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

    public static boolean areBasesEquivalent(List<Polynomial> setA, List<Polynomial> setB) {
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

        // arg0 = number of polynomials to generate
        // arg1 = term order: 0 = Lex, 1 = GrLex, 2 = RevLex

        if (mode != 0) {
            int numPolynomials = args.length > 0 ? Integer.parseInt(args[0]) : 3;
            LCG rand = new LCG(12345, 1345, 16645, 1013904);
            List<Polynomial> inputBasis = new ArrayList<>();
            GrobnerSmart.termOrder = TermOrder.Lex;
            if (args.length > 1) {
                int orderArg = Integer.parseInt(args[1]);
                switch (orderArg) {
                    case 0: GrobnerSmart.termOrder = TermOrder.Lex; break;
                    case 1: GrobnerSmart.termOrder = TermOrder.GrLex; break;
                    case 2: GrobnerSmart.termOrder = TermOrder.RevLex; break;
                    default: GrobnerSmart.termOrder = TermOrder.Lex;
                }
            }
            for (int i = 0; i < numPolynomials; i++) {
                int numTerms = 3;
                List<Term> terms = new ArrayList<>();
                for (int j = 0; j < numTerms; j++) {
                    long coefficient = rand.nextInt() % modulus;
                    int[] exps = new int[6];
                    for (int k = 0; k < 3; k++) {
                        exps[k] = rand.nextInt() % 4;
                    }
                    terms.add(Term.fromExponents(coefficient, exps));
                }
                Polynomial poly = new Polynomial(terms);
                inputBasis.add(poly);
            }
            List<Polynomial> basis = naiveGrobnerBasis(inputBasis);
            System.out.println("Computed GrobnerSmart Basis Polynomials:");
            for (Polynomial poly : basis) {
                System.out.println(poly);
            }
        }
        
        else {
            modulus = 7;
            GrobnerSmart.termOrder = TermOrder.Lex;
            int n = args.length > 0 ? Integer.parseInt(args[0]) : 6;
            if (n == 4) {
                System.out.println("Java Specialized finite coeff bitpacked exp cyclic 4");  
                // Cyclic 4 system: x0+x1+x2+x3, x0x1+x1x2+x2x3+x3x0, x0x1x2+x1x2x3+x2x3x0+x3x0x1, x0x1x2x3-1
                // f1 = x0 + x1 + x2 + x3
                List<Term> terms1 = new ArrayList<>();
                terms1.add(Term.fromExponents(1L, new int[]{1, 0, 0, 0, 0, 0}));
                terms1.add(Term.fromExponents(1L, new int[]{0, 1, 0, 0, 0, 0}));
                terms1.add(Term.fromExponents(1L, new int[]{0, 0, 1, 0, 0, 0}));
                terms1.add(Term.fromExponents(1L, new int[]{0, 0, 0, 1, 0, 0}));
                Polynomial p1 = new Polynomial(terms1);

                // f2 = x0x1 + x1x2 + x2x3 + x3x0
                List<Term> terms2 = new ArrayList<>();
                terms2.add(Term.fromExponents(1L, new int[]{1, 1, 0, 0, 0, 0}));
                terms2.add(Term.fromExponents(1L, new int[]{0, 1, 1, 0, 0, 0}));
                terms2.add(Term.fromExponents(1L, new int[]{0, 0, 1, 1, 0, 0}));
                terms2.add(Term.fromExponents(1L, new int[]{1, 0, 0, 1, 0, 0}));
                Polynomial p2 = new Polynomial(terms2);

                // f3 = x0x1x2 + x1x2x3 + x2x3x0 + x3x0x1
                List<Term> terms3 = new ArrayList<>();
                terms3.add(Term.fromExponents(1L, new int[]{1, 1, 1, 0, 0, 0}));
                terms3.add(Term.fromExponents(1L, new int[]{0, 1, 1, 1, 0, 0}));
                terms3.add(Term.fromExponents(1L, new int[]{1, 0, 1, 1, 0, 0}));
                terms3.add(Term.fromExponents(1L, new int[]{1, 1, 0, 1, 0, 0}));
                Polynomial p3 = new Polynomial(terms3);

                // f4 = x0x1x2x3 - 1
                List<Term> terms4 = new ArrayList<>();
                terms4.add(Term.fromExponents(1L, new int[]{1, 1, 1, 1, 0, 0}));
                terms4.add(Term.fromExponents(modulus-1L, new int[]{0, 0, 0, 0, 0, 0}));
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
            else if (n == 5 ) {
                System.out.println("Java Specialized finite coeff bitpacked exp cyclic 5");
                List<Term> terms1 = new ArrayList<>();
                terms1.add(Term.fromExponents(1L, new int[]{1, 0, 0, 0, 0, 0}));
                terms1.add(Term.fromExponents(1L, new int[]{0, 1, 0, 0, 0, 0}));
                terms1.add(Term.fromExponents(1L, new int[]{0, 0, 1, 0, 0, 0}));
                terms1.add(Term.fromExponents(1L, new int[]{0, 0, 0, 1, 0, 0}));
                terms1.add(Term.fromExponents(1L, new int[]{0, 0, 0, 0, 1, 0}));
                Polynomial p1 = new Polynomial(terms1);

                List<Term> terms2 = new ArrayList<>();
                terms2.add(Term.fromExponents(1L, new int[]{1, 1, 0, 0, 0, 0}));
                terms2.add(Term.fromExponents(1L, new int[]{0, 1, 1, 0, 0, 0}));
                terms2.add(Term.fromExponents(1L, new int[]{0, 0, 1, 1, 0, 0}));
                terms2.add(Term.fromExponents(1L, new int[]{0, 0, 0, 1, 1, 0}));
                terms2.add(Term.fromExponents(1L, new int[]{1, 0, 0, 0, 1, 0}));
                Polynomial p2 = new Polynomial(terms2);

                List<Term> terms3 = new ArrayList<>();
                terms3.add(Term.fromExponents(1L, new int[]{1, 1, 1, 0, 0, 0}));
                terms3.add(Term.fromExponents(1L, new int[]{0, 1, 1, 1, 0, 0}));
                terms3.add(Term.fromExponents(1L, new int[]{0, 0, 1, 1, 1, 0}));
                terms3.add(Term.fromExponents(1L, new int[]{1, 0, 0, 1, 1, 0}));
                terms3.add(Term.fromExponents(1L, new int[]{1, 1, 0, 0, 1, 0}));
                Polynomial p3 = new Polynomial(terms3); 

                List<Term> terms4 = new ArrayList<>();
                terms4.add(Term.fromExponents(1L, new int[]{1, 1, 1, 1, 0, 0}));
                terms4.add(Term.fromExponents(1L, new int[]{0, 1, 1, 1, 1, 0}));
                terms4.add(Term.fromExponents(1L, new int[]{1, 0, 1, 1, 1, 0}));
                terms4.add(Term.fromExponents(1L, new int[]{1, 1, 0, 1, 1, 0}));
                terms4.add(Term.fromExponents(1L, new int[]{1, 1, 1, 0, 1, 0}));
                Polynomial p4 = new Polynomial(terms4);

                List<Term> terms5 = new ArrayList<>();
                terms5.add(Term.fromExponents(1L, new int[]{1, 1, 1, 1, 1, 0}));
                terms5.add(Term.fromExponents(modulus-1L, new int[]{0, 0, 0, 0, 0, 0}));
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
                System.out.println("Java Specialized finite coeff bitpacked exp cyclic 6");
                List<Term> terms1 = new ArrayList<>();
                terms1.add(Term.fromExponents(1L, new int[]{1, 0, 0, 0, 0, 0}));
                terms1.add(Term.fromExponents(1L, new int[]{0, 1, 0, 0, 0, 0}));
                terms1.add(Term.fromExponents(1L, new int[]{0, 0, 1, 0, 0, 0}));
                terms1.add(Term.fromExponents(1L, new int[]{0, 0, 0, 1, 0, 0}));
                terms1.add(Term.fromExponents(1L, new int[]{0, 0, 0, 0, 1, 0}));
                terms1.add(Term.fromExponents(1L, new int[]{0, 0, 0, 0, 0, 1}));
                Polynomial p1 = new Polynomial(terms1);

                List<Term> terms2 = new ArrayList<>();
                terms2.add(Term.fromExponents(1L, new int[]{1, 1, 0, 0, 0, 0}));
                terms2.add(Term.fromExponents(1L, new int[]{0, 1, 1, 0, 0, 0}));
                terms2.add(Term.fromExponents(1L, new int[]{0, 0, 1, 1, 0, 0}));
                terms2.add(Term.fromExponents(1L, new int[]{0, 0, 0, 1, 1, 0}));
                terms2.add(Term.fromExponents(1L, new int[]{0, 0, 0, 0, 1, 1}));
                terms2.add(Term.fromExponents(1L, new int[]{1, 0, 0, 0, 0, 1}));
                Polynomial p2 = new Polynomial(terms2);

                List<Term> terms3 = new ArrayList<>();
                terms3.add(Term.fromExponents(1L, new int[]{1, 1, 1, 0, 0, 0}));
                terms3.add(Term.fromExponents(1L, new int[]{0, 1, 1, 1, 0, 0}));
                terms3.add(Term.fromExponents(1L, new int[]{0, 0, 1, 1, 1, 0}));
                terms3.add(Term.fromExponents(1L, new int[]{0, 0, 0, 1, 1, 1}));
                terms3.add(Term.fromExponents(1L, new int[]{1, 0, 0, 0, 1, 1}));
                terms3.add(Term.fromExponents(1L, new int[]{1, 1, 0, 0, 0, 1}));
                Polynomial p3 = new Polynomial(terms3);

                List<Term> terms4 = new ArrayList<>();
                terms4.add(Term.fromExponents(1L, new int[]{1, 1, 1, 1, 0, 0}));
                terms4.add(Term.fromExponents(1L, new int[]{0, 1, 1, 1, 1, 0}));
                terms4.add(Term.fromExponents(1L, new int[]{0, 0, 1, 1, 1, 1}));
                terms4.add(Term.fromExponents(1L, new int[]{1, 0, 0, 1, 1, 1}));
                terms4.add(Term.fromExponents(1L, new int[]{1, 1, 0, 0, 1, 1}));
                terms4.add(Term.fromExponents(1L, new int[]{1, 1, 1, 0, 0, 1}));
                Polynomial p4 = new Polynomial(terms4);

                List<Term> terms5 = new ArrayList<>();
                terms5.add(Term.fromExponents(1L, new int[]{1, 1, 1, 1, 1, 0}));
                terms5.add(Term.fromExponents(1L, new int[]{0, 1, 1, 1, 1, 1}));
                terms5.add(Term.fromExponents(1L, new int[]{1, 0, 1, 1, 1, 1}));
                terms5.add(Term.fromExponents(1L, new int[]{1, 1, 0, 1, 1, 1}));
                terms5.add(Term.fromExponents(1L, new int[]{1, 1, 1, 0, 1, 1}));
                terms5.add(Term.fromExponents(1L, new int[]{1, 1, 1, 1, 0, 1}));
                Polynomial p5 = new Polynomial(terms5);

                List<Term> terms6 = new ArrayList<>();
                terms6.add(Term.fromExponents(1L, new int[]{1, 1, 1, 1, 1, 1}));
                terms6.add(Term.fromExponents(modulus-1L, new int[]{0, 0, 0, 0, 0, 0}));
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

    // Helper static methods for bitpacking
    public static int[] unpackExponents(long packed) {
        int[] exps = new int[6];
        for (int i = 0; i < 6; i++) {
            exps[i] = (int) ((packed >> (40 - 8 * i)) & 0xFF);
        }
        return exps;
    }
    public static long packExponents(int[] exps) {
        long packed = 0L;
        int degree = 0;
        for (int i = 0; i < 6; i++) {
            packed |= ((long) (exps[i] & 0xFF)) << (40 - 8 * i);
            degree += exps[i];
        }
        packed |= ((long) degree) << 48;
        return packed;
    }
}
