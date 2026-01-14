
package generic;

import java.util.*;
import helpers.LCG;


public class GenGrobner {
    public enum TermOrder {
        Lex, GrLex, RevLex
    }

    public static class Term<C extends IField<C>, E extends IExponent<E>> {
        public final C coefficient;
        public final E exponents;

        public Term(C coefficient, E exponents) {
            this.coefficient = coefficient;
            this.exponents = exponents;
        }

        public int degree() {
            return exponents.degree();
        }

        public int compareTo(Term<C, E> other, TermOrder order) {
            switch (order) {
                case Lex:
                    return exponents.compareTo(other.exponents);
                case GrLex:
                    int selfDegree = this.degree();
                    int otherDegree = other.degree();
                    int cmp = Integer.compare(selfDegree, otherDegree);
                    if (cmp != 0) return cmp;
                    return exponents.compareTo(other.exponents);
                case RevLex:
                    selfDegree = this.degree();
                    otherDegree = other.degree();
                    cmp = Integer.compare(selfDegree, otherDegree);
                    if (cmp != 0) return cmp;
                    return other.exponents.compareTo(exponents);
                default:
                    return exponents.compareTo(other.exponents);
            }
        }

        public boolean canReduce(Term<C, E> divisorLeading) {
            return exponents.canReduce(divisorLeading.exponents);
        }

        public E lcm(Term<C, E> other) {
            return exponents.lcm(other.exponents);
        }

        @Override
        public boolean equals(Object o) {
            if (this == o) return true;
            if (!(o instanceof Term)) return false;
            Term<?, ?> term = (Term<?, ?>) o;
            return Objects.equals(coefficient, term.coefficient) && Objects.equals(exponents, term.exponents);
        }

        @Override
        public int hashCode() {
            return Objects.hash(coefficient, exponents);
        }

        @Override
        public String toString() {
            return coefficient + "*" + exponents;
        }
    }

    public static class Polynomial<C extends IField<C>, E extends IExponent<E>> {
        public final List<Term<C, E>> terms;

        public Polynomial(List<Term<C, E>> terms, TermOrder order) {
            List<Term<C, E>> filtered = new ArrayList<>();
            for (Term<C, E> t : terms) {
                if (Math.abs(t.coefficient.coerce() - 0.0) > 1e-2) filtered.add(t);
            }
            filtered.sort((a, b) -> b.compareTo(a, order));
            this.terms = Collections.unmodifiableList(filtered);
        }

        public Polynomial<C,E> deepCopy(TermOrder order) {
            List<Term<C, E>> copiedTerms = new ArrayList<>();
            for (Term<C, E> t : terms) {
                copiedTerms.add(new Term<>(t.coefficient, t.exponents));
            }
            return new Polynomial<>(copiedTerms, order);
        }

        public Polynomial<C, E> add(Polynomial<C, E> other, TermOrder order) {
            //List<Term<C, E>> result = this.deepCopy(order).terms;
            List<Term<C, E>> result = new ArrayList<>(terms);
            for (Term<C, E> t : other.terms) {
                boolean found = false;
                for (int i = 0; i < result.size(); i++) {
                    Term<C, E> rt = result.get(i);
                    if (rt.exponents.equals(t.exponents)) {
                        // making a new term to avoid mutating existing ones
                        result.set(i, new Term<>(rt.coefficient.a(t.coefficient), rt.exponents));
                        found = true;
                        break;
                    }
                }
                if (!found) result.add(t);
            }
            return new Polynomial<>(result, order);
        }

        public Polynomial<C, E> subtract(Polynomial<C, E> other, TermOrder order) {
            List<Term<C, E>> result = new ArrayList<>(terms);
            for (Term<C, E> t : other.terms) {
                boolean found = false;
                for (int i = 0; i < result.size(); i++) {
                    Term<C, E> rt = result.get(i);
                    if (rt.exponents.equals(t.exponents)) {
                        // making a new term to avoid mutating existing ones
                        result.set(i, new Term<>(rt.coefficient.s(t.coefficient), rt.exponents));
                        found = true;
                        break;
                    }
                }
                if (!found) result.add(new Term<>(t.coefficient.zero().s(t.coefficient), t.exponents));
            }
            return new Polynomial<>(result, order);
        }

        public Polynomial<C, E> multiplyByTerm(Term<C, E> term, TermOrder order) {
            List<Term<C, E>> result = new ArrayList<>();
            for (Term<C, E> t : terms) {
                result.add(new Term<>(t.coefficient.m(term.coefficient), t.exponents.add(term.exponents)));
            }
            return new Polynomial<>(result, order);
        }

        public Polynomial<C, E> reduce(List<Polynomial<C, E>> divisors, TermOrder order) {
            Polynomial<C, E> result = this;
            boolean reduced;
            do {
                reduced = false;
                for (Polynomial<C, E> divisor : divisors) {
                    if (result.terms.isEmpty() || divisor.terms.isEmpty()) continue;
                    Term<C, E> lead = result.terms.get(0);
                    Term<C, E> divLead = divisor.terms.get(0);
                    //System.out.println("Attempting to reduce " + lead + " by " + divLead);
                    if (lead.canReduce(divLead)) {
                        C coeff = lead.coefficient.d(divLead.coefficient);
                        E exps = lead.exponents.sub(divLead.exponents);
                        Term<C, E> reductionTerm = new Term<>(coeff, exps);
                        Polynomial<C, E> scaledDivisor = divisor.multiplyByTerm(reductionTerm, order);
                        result = result.subtract(scaledDivisor, order);
                        reduced = true;
                        //System.out.println("Reduced to: " + result + "\n\n");
                        break;
                    }
                }
            } while (reduced);
            return new Polynomial<>(result.terms, order);
        }

        public static <C extends IField<C>, E extends IExponent<E>> Polynomial<C, E> sPolynomial(Polynomial<C, E> p1, Polynomial<C, E> p2, TermOrder order) {
            Term<C, E> lead1 = p1.terms.get(0);
            Term<C, E> lead2 = p2.terms.get(0);
            E lcmExps = lead1.lcm(lead2);
            E scale1 = lcmExps.sub(lead1.exponents);
            E scale2 = lcmExps.sub(lead2.exponents);
            Polynomial<C, E> scaled1 = p1.multiplyByTerm(new Term<>(lead1.coefficient.one(), scale1), order);
            Polynomial<C, E> scaled2 = p2.multiplyByTerm(new Term<>(lead2.coefficient.one(), scale2), order);
            return scaled1.subtract(scaled2, order);
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
    }

    public static <C extends IField<C>, E extends IExponent<E>> List<Polynomial<C, E>> naiveGrobnerBasis(List<Polynomial<C, E>> polys, TermOrder order) {
        List<Polynomial<C, E>> basis = new ArrayList<>(polys);
        Set<Polynomial<C, E>> basisSet = new HashSet<>(basis);
        boolean added;
        do {
            added = false;
            int n = basis.size();
            for (int i = 0; i < n; i++) {
                for (int j = i + 1; j < n; j++) {
                    Polynomial<C, E> sPoly = Polynomial.sPolynomial(basis.get(i), basis.get(j), order);
                    Polynomial<C, E> reduced = sPoly.reduce(basis, order);
                    if (!reduced.terms.isEmpty() && !basisSet.contains(reduced)) {
                        //System.out.println("Reduced S-Polynomial of basis[" + i + "] and basis[" + j + "]: " + reduced);
                        basisSet.add(reduced);
                        basis.add(reduced);
                        added = true;
                    }
                }
            }
        } while (added);
        // print basis before reduction
        //System.out.println("Basis before reduction:");
        /*for (Polynomial<C, E> poly : basis) {
            System.out.println(poly);
        }*/

        List<Polynomial<C, E>> reducedBasis = new ArrayList<>();
        for (Polynomial<C, E> poly : basis) {
            List<Polynomial<C, E>> basisExcludingSelf = new ArrayList<>(basis);
            basisExcludingSelf.remove(poly);
            Polynomial<C, E> reduced = poly.reduce(basisExcludingSelf, order);
            //System.out.println("Reducing");
            if (!reduced.terms.isEmpty() && !reducedBasis.contains(reduced)) {
                reducedBasis.add(reduced);
            }
        }
        return reducedBasis;
    }

    // Example usage for SingleField + VecExponent
    public static void main(String[] args) {
        int mode = 0;
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
        // print ordering
        System.out.println("Using term order: " + order);
        if (mode != 0) {
            int numPolynomials = args.length > 0 ? Integer.parseInt(args[0]) : 3;
            LCG rand = new LCG(12345, 1345, 65, 17);
            List<Polynomial<SingleField, VecExponent>> inputBasis = new ArrayList<>();
            for (int i = 0; i < numPolynomials; i++) {
                List<Term<SingleField, VecExponent>> terms = new ArrayList<>();
                for (int j = 0; j < 3; j++) {
                    SingleField coefficient = new SingleField((float)(rand.nextDouble() * 2.0 - 1.0));
                    VecExponent exponents = new VecExponent(Arrays.asList(rand.nextInt() % 4, rand.nextInt() % 4, rand.nextInt() % 4));
                    terms.add(new Term<>(coefficient, exponents));
                }
                inputBasis.add(new Polynomial<>(terms, order));
            }
            List<Polynomial<SingleField, VecExponent>> basis = naiveGrobnerBasis(inputBasis, order);
            System.out.println(basis.size());
            System.out.println("Computed Grobner Basis Polynomials:");
            for (Polynomial<SingleField, VecExponent> poly : basis) {
                System.out.println(poly);
            }
        } else {
            // Example: x^2*y + y^2*z + z^2*x
            List<Term<SingleField, BitPackedExponent>> terms1 = new ArrayList<>();
            terms1.add(new Term<>(new SingleField(1.0f), BitPackedExponent.fromArray(new int[]{2,1,0,0,0,0})));
            terms1.add(new Term<>(new SingleField(1.0f), BitPackedExponent.fromArray(new int[]{0,2,1,0,0,0})));
            terms1.add(new Term<>(new SingleField(1.0f), BitPackedExponent.fromArray(new int[]{1,0,2,0,0,0})));
            Polynomial<SingleField, BitPackedExponent> p1 = new Polynomial<>(terms1, order);
            // xyz - 1
            List<Term<SingleField, BitPackedExponent>> terms2 = new ArrayList<>();
            terms2.add(new Term<>(new SingleField(1.0f), BitPackedExponent.fromArray(new int[]{1,1,1,0,0,0})));
            terms2.add(new Term<>(new SingleField(-1.0f), BitPackedExponent.fromArray(new int[]{0,0,0,0,0,0})));
            Polynomial<SingleField, BitPackedExponent> p2 = new Polynomial<>(terms2, order);
            // x + y + z
            List<Term<SingleField, BitPackedExponent>> terms3 = new ArrayList<>();
            terms3.add(new Term<>(new SingleField(1.0f), BitPackedExponent.fromArray(new int[]{1,0,0,0,0,0})));
            terms3.add(new Term<>(new SingleField(1.0f), BitPackedExponent.fromArray(new int[]{0,1,0,0,0,0})));
            terms3.add(new Term<>(new SingleField(1.0f), BitPackedExponent.fromArray(new int[]{0,0,1,0,0,0})));
            Polynomial<SingleField, BitPackedExponent> p3 = new Polynomial<>(terms3, order);
            List<Polynomial<SingleField, BitPackedExponent>> inputBasis = Arrays.asList(p1, p2, p3);
            System.out.println("Begin");


            // PRINT P1
            System.out.println("Input Polynomials:");
            for (Polynomial<SingleField, BitPackedExponent> poly : inputBasis) {
                System.out.println(poly);
            }
            System.out.println("Spoly of p1-p2: " + Polynomial.sPolynomial(p1, p2, order));

            List<Polynomial<SingleField, BitPackedExponent>> basis = naiveGrobnerBasis(inputBasis, order);
            System.out.println("Computed Grobner Basis Polynomials:");
            for (Polynomial<SingleField, BitPackedExponent> poly : basis) {
                System.out.println("\n\nPolynomial: " + poly);
            }
        }
    }
}
