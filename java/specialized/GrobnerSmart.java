package specialized;

import java.util.*;

public class GrobnerSmart {
    public static class Term {
        public double coefficient;
        public long exponents; // Bitpacked: [63..48]=degree, [47..40]=e0, [39..32]=e1, ... [7..0]=e5

        public Term(double coefficient, long exponents) {
            this.coefficient = coefficient;
            this.exponents = exponents;
        }

        // Create a term from array of exponents (length 6)
        public static Term fromExponents(double coefficient, int[] exps) {
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
        public int compareTo(Term other, TermOrder order) {
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
            return Math.abs(coefficient - term.coefficient) < 1e-5 && exponents == term.exponents;
        }

        @Override
        public int hashCode() {
            return Objects.hash(Math.round(coefficient * 1e5), exponents);
        }

        @Override
        public String toString() {
            return String.format("%.5f * %s", coefficient, Arrays.toString(unpackExponents(exponents)));
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

        public Polynomial(List<Term> terms, TermOrder order) {
            this.terms = new ArrayList<>(terms);
            this.terms.removeIf(t -> Math.abs(t.coefficient) < 1e-2);
            for (Term t : this.terms) {
                t.coefficient = Math.round(t.coefficient * 1e5) / 1e5;
            }
            //System.out.println("Before sort: " + this.terms);
            this.terms.sort((a, b) -> b.compareTo(a, order));
            //System.out.println("After sort: " + this.terms);
        }

        public Polynomial deepCopy(TermOrder order) {
            List<Term> newTerms = new ArrayList<>();
            for (Term t : this.terms) {
                newTerms.add(new Term(t.coefficient, t.exponents));
            }
            return new Polynomial(newTerms, order);
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

        public Polynomial add(Polynomial other, TermOrder order) {
            List<Term> result = new ArrayList<>(terms);
            for (Term term : other.terms) {
                boolean found = false;
                for (Term resTerm : result) {
                    if (resTerm.exponents == term.exponents) {
                        resTerm.coefficient += term.coefficient;
                        found = true;
                        break;
                    }
                }
                if (!found) {
                    result.add(new Term(term.coefficient, term.exponents));
                }
            }
            return new Polynomial(result, order);
        }

        public Polynomial subtract(Polynomial other, TermOrder order) {
            List<Term> result = new ArrayList<>(terms);
            for (Term term : other.terms) {
                boolean found = false;
                for (Term resTerm : result) {
                    if (resTerm.exponents == term.exponents) {
                        resTerm.coefficient -= term.coefficient;
                        found = true;
                        break;
                    }
                }
                if (!found) {
                    result.add(new Term(-term.coefficient, term.exponents));
                }
            }
            return new Polynomial(result, order);
        }

        public Polynomial multiplyByTerm(Term term, TermOrder order) {
            List<Term> result = new ArrayList<>();
            for (Term t : this.terms) {
                long exp = t.exponents + term.exponents;
                double coefficient = t.coefficient * term.coefficient;
                result.add(new Term(coefficient, exp));
            }
            return new Polynomial(result, order);
        }

        public Polynomial reduce(List<Polynomial> divisors, TermOrder order) {
            // deep copy of self
            Polynomial result = this.deepCopy(order);
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
                        double coefficient = leadingTerm.coefficient / divisorLeadingTerm.coefficient;
                        Term reductionTerm = new Term(coefficient, exponents);
                        Polynomial scaledDivisor = divisor.multiplyByTerm(reductionTerm, order);
                        result = result.subtract(scaledDivisor, order);
                        reduced = true;
                        //System.out.println("Reduced polynomial: " + result.toString());
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
            long lcmExponents = Term.lcm(leadingTermP1.exponents, leadingTermP2.exponents);
            long exp1 = lcmExponents - leadingTermP1.exponents;
            long exp2 = lcmExponents - leadingTermP2.exponents;
            Term multTerm1 = new Term(1.0, exp1);
            Term multTerm2 = new Term(1.0, exp2);
            Polynomial term1 = p1.multiplyByTerm(multTerm1, order);
            Polynomial term2 = p2.multiplyByTerm(multTerm2, order);
            return term1.subtract(term2, order);
            
        }
    }

    public static List<Polynomial> naiveGrobnerBasis(List<Polynomial> polynomials, TermOrder order) {
        // deep copies of input polynomials
        List<Polynomial> basis = new ArrayList<>();
        for (Polynomial poly : polynomials) {
            basis.add(poly.deepCopy(order));
        }
        Set<Polynomial> basisSet = new HashSet<>(basis);
        while (true) {
            //System.out.println("Grobner Basis Iteration, basis size: " + basis.size());
            boolean added = false;
            int basisLen = basis.size();
            for (int i = 0; i < basisLen; i++) {
                //System.out.println("Processing basis polynomial " + i);
                for (int j = i + 1; j < basisLen; j++) {
                    Polynomial sPoly = Polynomial.sPolynomial(basis.get(i), basis.get(j), order);
                    Polynomial reduced = sPoly.reduce(basis, order);
                    // print reduced
                    //System.out.println("Reduced S-Polynomial:" + reduced.toString());
                    if (!reduced.terms.isEmpty() && !basisSet.contains(reduced)) {
                        basis.add(reduced.deepCopy(order));
                        basisSet.add(reduced.deepCopy(order));
                        added = true;
                    }
                }
            }
           // System.out.println(added);
            if (!added) break;
        }
        //System.out.println("Checkpoint");
        List<Polynomial> reducedBasis = new ArrayList<>();
        
        for (Polynomial poly : basis) {
            //System.out.println("Basis: " + basis.toString());
            List<Polynomial> basisExcludingSelf = new ArrayList<>();
            for (Polynomial p : basis) {
                if (!p.equals(poly)) {
                    basisExcludingSelf.add(p.deepCopy(order));
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
        TermOrder order = TermOrder.Lex;
        // Example: x^2*y + y^2*z + z^2*x
        Term t1 = Term.fromExponents(1.0, new int[]{2, 1, 0, 0, 0, 0});
        Term t2 = Term.fromExponents(1.0, new int[]{0, 2, 1, 0, 0, 0});
        Term t3 = Term.fromExponents(1.0, new int[]{1, 0, 2, 0, 0, 0});
        Polynomial p1 = new Polynomial(Arrays.asList(t1, t2, t3), order);
        // x*y*z - 1
        Term t4 = Term.fromExponents(1.0, new int[]{1, 1, 1, 0, 0, 0});
        Term t5 = Term.fromExponents(-1.0, new int[]{0, 0, 0, 0, 0, 0});
        Polynomial p2 = new Polynomial(Arrays.asList(t4, t5), order);
        // x + y + z
        Term t6 = Term.fromExponents(1.0, new int[]{1, 0, 0, 0, 0, 0});
        Term t7 = Term.fromExponents(1.0, new int[]{0, 1, 0, 0, 0, 0});
        Term t8 = Term.fromExponents(1.0, new int[]{0, 0, 1, 0, 0, 0});
        Polynomial p3 = new Polynomial(Arrays.asList(t6, t7, t8), order);
        List<Polynomial> inputBasis = Arrays.asList(p1, p2, p3);
        List<Polynomial> basis = naiveGrobnerBasis(inputBasis, order);
        System.out.println("Computed GrobnerSmart Basis Polynomials:");
        for (Polynomial poly : basis) {
            System.out.println(poly);
        }

        // test removal of small coefficients
            Term t9 = Term.fromExponents(1.0, new int[]{2, 0, 0, 0, 0, 0});
            Term t10 = Term.fromExponents(1e-3, new int[]{1, 1, 0, 0, 0, 0});
            Polynomial p4 = new Polynomial(Arrays.asList(t9, t10), order);
            System.out.println("Polynomial with small coefficient term:");
            System.out.println(p4.toString());

        // Test S-poly
        Polynomial sPoly = Polynomial.sPolynomial(p1, p2, order);
        System.out.println("S-Polynomial of p1 and p2:");
        System.out.println(sPoly.toString());
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
