import { IField } from './iField';
import { IMath } from './iMath';
import { IOrdered } from './iOrdered';
import { ICopiable } from './iCopiable';
import { Term, Polynomial } from './polynomial';

export enum TermOrder {
  Lex,
  GrLex,
  RevLex
}

export function compareExponents<E>(a: E, b: E, order: TermOrder): number {
  // E must have deg() and lexCompare()
  // @ts-ignore
  switch (order) {
    case TermOrder.Lex:
      return a.lexCompare(b);
    case TermOrder.GrLex:
      // @ts-ignore
      const tDeg = a.deg();
      // @ts-ignore
      const oDeg = b.deg();
      if (tDeg !== oDeg) return tDeg - oDeg;
      return a.lexCompare(b);
    case TermOrder.RevLex:
      // @ts-ignore
      const tDeg2 = a.deg();
      // @ts-ignore
      const oDeg2 = b.deg();
      if (tDeg2 !== oDeg2) return tDeg2 - oDeg2;
      return -a.lexCompare(b);
    default:
      return 0;
  }
}

export function newPolynomial<C, E>(terms: Term<C, E>[], order: TermOrder): Polynomial<C, E> {
  // Remove near-zero coefficients and sort
  const filtered = terms.filter(t => Math.abs((t.coefficient as any).coerceToFloat()) > 1e-2);
  filtered.sort((a, b) => compareExponents(a.exponents, b.exponents, order));
  return new Polynomial(filtered);
}

export function addPolynomials<C, E>(a: Polynomial<C, E>, b: Polynomial<C, E>, order: TermOrder): Polynomial<C, E> {
  const result = [...a.terms];
  for (const t of b.terms) {
    let found = false;
    for (let i = 0; i < result.length; i++) {
      // @ts-ignore
      if (result[i].exponents.equals(t.exponents)) {
        // @ts-ignore
        result[i].coefficient = (result[i].coefficient as any).add(t.coefficient);
        found = true;
        break;
      }
    }
    if (!found) result.push(t);
  }
  return newPolynomial(result, order);
}

export function subtractPolynomials<C, E>(a: Polynomial<C, E>, b: Polynomial<C, E>, order: TermOrder): Polynomial<C, E> {
  const result = [...a.terms];
  for (const t of b.terms) {
    let found = false;
    for (let i = 0; i < result.length; i++) {
      // @ts-ignore
      if (result[i].exponents.equals(t.exponents)) {
        // @ts-ignore
        result[i].coefficient = (result[i].coefficient as any).sub(t.coefficient);
        found = true;
        break;
      }
    }
    if (!found) {
      // @ts-ignore
      result.push(new Term((t.coefficient as any).zero().sub(t.coefficient), t.exponents));
    }
  }
  return newPolynomial(result, order);
}

export function multiplyByTerm<C, E>(p: Polynomial<C, E>, term: Term<C, E>, order: TermOrder): Polynomial<C, E> {
  const terms = p.terms.map(t => new Term(
    // @ts-ignore
    (t.coefficient as any).mul(term.coefficient),
    // @ts-ignore
    t.exponents.add(term.exponents)
  ));
  return newPolynomial(terms, order);
}

export function reducePolynomial<C, E>(p: Polynomial<C, E>, divisors: Polynomial<C, E>[], order: TermOrder): Polynomial<C, E> {
  let result = p;
  while (true) {
    let reduced = false;
    if (result.terms.length === 0) break;
    const lead = result.terms[0];
    for (const divisor of divisors) {
      if (divisor.terms.length === 0) continue;
      const divLead = divisor.terms[0];
      // @ts-ignore
      if (lead.exponents.canReduce(divLead.exponents)) {
        // @ts-ignore
        const coeff = (lead.coefficient as any).div(divLead.coefficient);
        // @ts-ignore
        const exps = lead.exponents.sub(divLead.exponents);
        const reductionTerm = new Term(coeff, exps);
        const scaledDivisor = multiplyByTerm(divisor, reductionTerm, order);
        result = subtractPolynomials(result, scaledDivisor, order);
        reduced = true;
        break;
      }
    }
    if (!reduced) break;
  }
  return newPolynomial(result.terms, order);
}

export function sPolynomial<C, E>(p1: Polynomial<C, E>, p2: Polynomial<C, E>, order: TermOrder): Polynomial<C, E> {
  const lead1 = p1.terms[0];
  const lead2 = p2.terms[0];
  // @ts-ignore
  const lcmExps = lead1.exponents.lcm(lead2.exponents);
  // @ts-ignore
  const scale1 = lcmExps.sub(lead1.exponents);
  // @ts-ignore
  const scale2 = lcmExps.sub(lead2.exponents);
  // @ts-ignore
  const scaled1 = multiplyByTerm(p1, new Term((lead1.coefficient as any).one(), scale1), order);
  // @ts-ignore
  const scaled2 = multiplyByTerm(p2, new Term((lead2.coefficient as any).one(), scale2), order);
  return subtractPolynomials(scaled1, scaled2, order);
}

export function naiveGrobnerBasis<C, E>(polys: Polynomial<C, E>[], order: TermOrder): Polynomial<C, E>[] {
  let basis = [...polys];
  const basisSet = new Set<string>();
  let added = true;
  while (added) {
    added = false;
    const n = basis.length;
    for (let i = 0; i < n; i++) {
      for (let j = i + 1; j < n; j++) {
        const sPoly = sPolynomial(basis[i], basis[j], order);
        const reduced = reducePolynomial(sPoly, basis, order);
        const key = JSON.stringify(reduced.terms);
        if (reduced.terms.length > 0 && !basisSet.has(key)) {
          basisSet.add(key);
          basis.push(reduced);
          added = true;
        }
      }
    }
  }
  // Reduce basis by self
  const reducedBasis: Polynomial<C, E>[] = [];
  for (const poly of basis) {
    const basisExcl = basis.filter(p => p !== poly);
    const reduced = reducePolynomial(poly, basisExcl, order);
    const key = JSON.stringify(reduced.terms);
    if (reduced.terms.length > 0 && !reducedBasis.some(rb => JSON.stringify(rb.terms) === key)) {
      reducedBasis.push(reduced);
    }
  }
  return reducedBasis;
}
