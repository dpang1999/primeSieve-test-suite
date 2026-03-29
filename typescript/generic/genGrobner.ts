import { IField } from "./iField";
import { IExponent } from "./iExponent";
import { ICopiable } from "./iCopiable";
import { IntModP } from './intModP.js';
import { VecExponents } from './vecExponents.js';
import { BitPackedExponents } from "./bitPackedExponents.js";
import { parse } from "node:path";

export type Term<C extends IField<C>, E extends IExponent<E>> = {
  coefficient: C;
  exponents: E;
};

export enum TermOrder {
  Lex,
  GrLex,
  RevLex,
}
let termOrder = TermOrder.Lex;


function compare<C extends IField<C>, E extends IExponent<E>>(term: Term<C, E>, other: Term<C, E>): number {
  switch (termOrder) {
    case TermOrder.Lex:
      return term.exponents.lexCompare(other.exponents);
    case TermOrder.GrLex: {
      const degA = term.exponents.deg();
      const degB = other.exponents.deg();
      if (degA !== degB) return degA - degB;
      return term.exponents.lexCompare(other.exponents);
    }
    case TermOrder.RevLex: {
      const degA = term.exponents.deg();
      const degB = other.exponents.deg();
      if (degA !== degB) return degA - degB;
      else { return other.exponents.lexCompare(term.exponents); }
    }
    default:
      return 0;
  }
}

function polyKey(terms: Term<any, any>[]): string {
  let out = "";
  for (let i = 0; i < terms.length; i++) {
    if (i > 0) out += ";";
    out += terms[i].coefficient + "|";
    out += terms[i].exponents + ", ";
  }
  return out;
}

function can_reduce<C extends IField<C>, E extends IExponent<E>>(term: Term<C, E>, other: Term<C, E>): boolean {
  return term.exponents.canReduce(other.exponents);
}

function lcm<C extends IField<C>, E extends IExponent<E>>(term: Term<C, E>, other: Term<C, E>): E {
  return term.exponents.lcm(other.exponents);
}


class Polynomial<C extends IField<C> & ICopiable<C>, E extends IExponent<E>> {
  terms: Term<C, E>[];
  constructor(terms: Term<C, E>[]) {
    // Remove near-zero coefficients, round, and sort
    this.terms = terms
      .filter(t => Math.abs(t.coefficient.coerce_to_number()) > 0.0)
      .sort((a, b) => -compare(a, b));
  }

  add(other: Polynomial<C, E>): Polynomial<C, E> {
    const results = this.terms.map(t => ({ coefficient: t.coefficient.copy(), exponents: t.exponents }));
    for (const t of other.terms) {
      let found = false;
      for (const rt of results) {
        if (rt.exponents.equals(t.exponents)) {
          rt.coefficient = rt.coefficient.a(t.coefficient);
          found = true;
          break;
        }
      }
      if (!found) results.push({ coefficient: t.coefficient, exponents: t.exponents });
    }
    return new Polynomial<C, E>(results);
  }

  subtract(other: Polynomial<C, E>): Polynomial<C, E> {
    const results = this.terms.map(t => ({ coefficient: t.coefficient.copy(), exponents: t.exponents }));
    for (const t of other.terms) {
      let found = false;
      for (const rt of results) {
        if (rt.exponents.equals(t.exponents)) {
          rt.coefficient = rt.coefficient.s(t.coefficient);
          found = true;
          break;
        }
      }
      if (!found) results.push({ coefficient: t.coefficient.zero().s(t.coefficient), exponents: t.exponents });
    }
    return new Polynomial(results);
  }

  multiplyByTerm(term: Term<C, E>): Polynomial<C, E> {
    const terms = this.terms.map(t => ({
      coefficient: t.coefficient.m(term.coefficient),
      exponents: t.exponents.add(term.exponents),
    }));
    return new Polynomial(terms);
  }

  reduce(divisors: Polynomial<C, E>[]): Polynomial<C, E> {
    let result = new Polynomial(this.terms);
    const remainder: Term<C, E>[] = [];
    while (true) {
      let reduced = false;
      for (const divisor of divisors) {
        if (result.terms.length === 0 || divisor.terms.length === 0) continue;
        const lead = result.terms[0];
        const divLead = divisor.terms[0];
        if (can_reduce(lead, divLead)) {
          const coeff = lead.coefficient.d(divLead.coefficient);
          const exps = lead.exponents.sub(divLead.exponents);
          const reductionTerm: Term<C,E> = { coefficient: coeff, exponents: exps };
          const scaledDivisor = divisor.multiplyByTerm(reductionTerm);
          result = result.subtract(scaledDivisor);
          reduced = true;
          break;
        }
      }
      if (!reduced) {
        if (result.terms.length === 0) break;
        remainder.push(result.terms[0]);
        result.terms = result.terms.slice(1);
      }
    }
    result.terms.push(...remainder);
    return new Polynomial(result.terms);
  }

   makeMonic(): Polynomial<C,E> {
      if (this.terms.length === 0) return this;
      const leadCoeff = this.terms[0].coefficient;
      if (leadCoeff.coerce_to_number() === 0) return this;
      const newTerms = this.terms.map(t => ({
        coefficient: (t.coefficient.d(leadCoeff)),
        exponents: t.exponents,
      }));
      return new Polynomial(newTerms);
    }

  static sPolynomial<C extends IField<C> & ICopiable<C>, E extends IExponent<E>>(p1: Polynomial<C, E>, p2: Polynomial<C, E>): Polynomial<C, E> {
    const lead1 = p1.terms[0];
    const lead2 = p2.terms[0];
    const lcmExps = lcm(lead1, lead2);
    const scale1 = lcmExps.sub(lead1.exponents);
    const scale2 = lcmExps.sub(lead2.exponents);
    const scaled1 = p1.multiplyByTerm({ coefficient: lead2.coefficient.copy(), exponents: scale1 });
    const scaled2 = p2.multiplyByTerm({ coefficient: lead1.coefficient.copy(), exponents: scale2 });
    return scaled1.subtract(scaled2);
  }

toString() {
    return this.terms.map(term => term.coefficient.toString() + "*" + term.exponents.toString()).join(" + ");
  }
}

export function naiveGrobnerBasis<C extends IField<C> & ICopiable<C>,  E extends IExponent<E>>(polys: Polynomial<C, E>[]): Polynomial<C, E>[] {
  let basis = polys.map(p => new Polynomial(p.terms));
  const basisSet = new Set<string>();

  for (const poly of basis) {
    basisSet.add(polyKey(poly.terms));
  }
  
  const pairs: Array<[number, number]> = [];
  for (let i = 0; i < basis.length; i++) {
    for (let j = i + 1; j < basis.length; j++) {
      pairs.push([i, j]);
    }
  }

  while (pairs.length > 0) {
    const [i, j] = pairs.shift()!;
    const sPoly = Polynomial.sPolynomial(basis[i], basis[j]);
    const reduced = sPoly.reduce(basis);
    const key = polyKey(reduced.terms);
    if (reduced.terms.length > 0 && !basisSet.has(key)) {
      basisSet.add(key);
      const newIdx = basis.length;
      basis.push(reduced);
      // Add new pairs with the new basis element
      for (let k = 0; k < newIdx; k++) {
        pairs.push([k, newIdx]);
      }
    }
  }
  
  // Reduce basis by self
  const reducedBasis: Polynomial<C, E>[] = [];
  for (const poly of basis) {
    const basisExcl = basis.filter(p => p !== poly);
    const reduced = poly.reduce(basisExcl);
    const key = polyKey(reduced.terms);
    if (reduced.terms.length > 0 && !reducedBasis.some(rb => polyKey(rb.terms) === key)) {
      reducedBasis.push(reduced.makeMonic());
    }
  }
  return reducedBasis;
}




function main() {
  const args = process.argv.slice(2);
  const mode = 0;
  if (mode != 0) {

  }
  else {
    // Test cyclic-4, cyclic-5, cyclic-6
    const n = args[0] ? parseInt(args[0]) : 4;
    const vecType = args[1] ? parseInt(args[1]) : 0;
    // 0 for vec exponent, 1 for bitpacked
    const modulus = 7;
    IntModP.setModulus(modulus);

    if (n == 4) {
      if (vecType == 0) {
        console.log("TypeScript generic finite coeff vec exponent cyclic 4");
        // a + b + c + d
        const q1 = new Polynomial<IntModP, VecExponents>([
          { coefficient: new IntModP(1), exponents: new VecExponents([1, 0, 0, 0]) },
          { coefficient: new IntModP(1), exponents: new VecExponents([0, 1, 0, 0]) },
          { coefficient: new IntModP(1), exponents: new VecExponents([0, 0, 1, 0]) },
          { coefficient: new IntModP(1), exponents: new VecExponents([0, 0, 0, 1]) },
        ]);
        // ab + bc + cd + ad
        const q2 = new Polynomial<IntModP, VecExponents>([
          { coefficient: new IntModP(1), exponents: new VecExponents([1, 1, 0, 0]) },
          { coefficient: new IntModP(1), exponents: new VecExponents([0, 1, 1, 0]) },
          { coefficient: new IntModP(1), exponents: new VecExponents([0, 0, 1, 1]) },
          { coefficient: new IntModP(1), exponents: new VecExponents([1, 0, 0, 1]) },
        ]);
        // abc + bcd + cda + dab
        const q3 = new Polynomial<IntModP, VecExponents>([
          { coefficient: new IntModP(1), exponents: new VecExponents([1, 1, 1, 0]) },
          { coefficient: new IntModP(1), exponents: new VecExponents([0, 1, 1, 1]) },
          { coefficient: new IntModP(1), exponents: new VecExponents([1, 0, 1, 1]) },
          { coefficient: new IntModP(1), exponents: new VecExponents([1, 1, 0, 1]) },
        ]);
        // abcd - 1
        const q4 = new Polynomial<IntModP, VecExponents>([
          { coefficient: new IntModP(1), exponents: new VecExponents([1, 1, 1, 1]) },
          { coefficient: new IntModP(modulus - 1), exponents: new VecExponents([0, 0, 0, 0]) },
        ]);
        const start = [q1, q2, q3, q4];
        for (let i = 0; i < 10; i++) {
          const basis = naiveGrobnerBasis(start);
          console.log("Iteration " + i + " complete");
            if (i == 9) {
            console.log(`Output basis size: ${basis.length}`);
            for (let i = 0; i < basis.length; i++) {
              console.log(`  [${i}]: ${basis[i].toString()}`);
            }
          }
        }
      }
      else {
        console.log("TypeScript generic finite coeff bitpacked exponent cyclic 4");
        // a + b + c + d
        const q1 = new Polynomial<IntModP, BitPackedExponents>([
          { coefficient: new IntModP(1), exponents: BitPackedExponents.fromArray([1, 0, 0, 0, 0, 0]) },
          { coefficient: new IntModP(1), exponents: BitPackedExponents.fromArray([0, 1, 0, 0, 0, 0]) },
          { coefficient: new IntModP(1), exponents: BitPackedExponents.fromArray([0, 0, 1, 0, 0, 0]) },
          { coefficient: new IntModP(1), exponents: BitPackedExponents.fromArray([0, 0, 0, 1, 0, 0]) },
        ]);
        // ab + bc + cd + ad
        const q2 = new Polynomial<IntModP, BitPackedExponents>([
          { coefficient: new IntModP(1), exponents: BitPackedExponents.fromArray([1, 1, 0, 0, 0, 0]) },
          { coefficient: new IntModP(1), exponents: BitPackedExponents.fromArray([0, 1, 1, 0, 0, 0]) },
          { coefficient: new IntModP(1), exponents: BitPackedExponents.fromArray([0, 0, 1, 1, 0, 0]) },
          { coefficient: new IntModP(1), exponents: BitPackedExponents.fromArray([1, 0, 0, 1, 0, 0]) },
        ]);
        // abc + bcd + cda + dab
        const q3 = new Polynomial<IntModP, BitPackedExponents>([
          { coefficient: new IntModP(1), exponents: BitPackedExponents.fromArray([1, 1, 1, 0, 0, 0]) },
          { coefficient: new IntModP(1), exponents: BitPackedExponents.fromArray([0, 1, 1, 1, 0, 0]) },
          { coefficient: new IntModP(1), exponents: BitPackedExponents.fromArray([1, 0, 1, 1, 0, 0]) },
          { coefficient: new IntModP(1), exponents: BitPackedExponents.fromArray([1, 1, 0, 1, 0, 0]) },
        ]);
        // abcd - 1
        const q4 = new Polynomial<IntModP, BitPackedExponents>([
          { coefficient: new IntModP(1), exponents: BitPackedExponents.fromArray([1, 1, 1, 1, 0, 0]) },
          { coefficient: new IntModP(modulus - 1), exponents: BitPackedExponents.fromArray([0, 0, 0, 0, 0, 0]) },
        ]);
        const start = [q1, q2, q3, q4];
        for (let i = 0; i < 10; i++) {
          const basis = naiveGrobnerBasis(start);
          console.log("Iteration " + i + " complete");
            if (i == 9) {
            console.log(`Output basis size: ${basis.length}`);
            for (let i = 0; i < basis.length; i++) {
              console.log(`  [${i}]: ${basis[i].toString()}`);
            }
          }
        }
      }
    }
    else if (n == 5){
      if (vecType == 0) {
        console.log("TypeScript generic finite coeff vec exponent cyclic 5");
        const p1 = new Polynomial<IntModP, VecExponents>([
          { coefficient: new IntModP(1), exponents: new VecExponents([1, 0, 0, 0, 0]) },
          { coefficient: new IntModP(1), exponents: new VecExponents([0, 1, 0, 0, 0]) },
          { coefficient: new IntModP(1), exponents: new VecExponents([0, 0, 1, 0, 0]) },
          { coefficient: new IntModP(1), exponents: new VecExponents([0, 0, 0, 1, 0]) },
          { coefficient: new IntModP(1), exponents: new VecExponents([0, 0, 0, 0, 1]) },
        ]);
        const p2 = new Polynomial<IntModP, VecExponents>([
          { coefficient: new IntModP(1), exponents: new VecExponents([1, 1, 0, 0, 0]) },
          { coefficient: new IntModP(1), exponents: new VecExponents([0, 1, 1, 0, 0]) },
          { coefficient: new IntModP(1), exponents: new VecExponents([0, 0, 1, 1, 0]) },
          { coefficient: new IntModP(1), exponents: new VecExponents([0, 0, 0, 1, 1]) },
          { coefficient: new IntModP(1), exponents: new VecExponents([1, 0, 0, 0, 1]) },
        ]);
        const p3 = new Polynomial<IntModP, VecExponents>([
          { coefficient: new IntModP(1), exponents: new VecExponents([1, 1, 1, 0, 0]) },
          { coefficient: new IntModP(1), exponents: new VecExponents([0, 1, 1, 1, 0]) },
          { coefficient: new IntModP(1), exponents: new VecExponents([0, 0, 1, 1, 1]) },
          { coefficient: new IntModP(1), exponents: new VecExponents([1, 0, 0, 1, 1]) },
          { coefficient: new IntModP(1), exponents: new VecExponents([1, 1, 0, 0, 1]) },
        ]);
        const p4 = new Polynomial<IntModP, VecExponents>([
          { coefficient: new IntModP(1), exponents: new VecExponents([1, 1, 1, 1, 0]) },
          { coefficient: new IntModP(1), exponents: new VecExponents([0, 1, 1, 1, 1]) },
          { coefficient: new IntModP(1), exponents: new VecExponents([1, 0, 1, 1, 1]) },
          { coefficient: new IntModP(1), exponents: new VecExponents([1, 1, 0, 1, 1]) },
          { coefficient: new IntModP(1), exponents: new VecExponents([1, 1, 1, 0, 1]) },
        ]);
        const p5 = new Polynomial<IntModP, VecExponents>([
          { coefficient: new IntModP(1), exponents: new VecExponents([1, 1, 1, 1, 1]) },
          { coefficient: new IntModP(modulus - 1), exponents: new VecExponents([0, 0, 0, 0, 0]) },
        ]);
        const start = [p1, p2, p3, p4, p5];
        for (let i = 0; i < 10; i++) {
          const basis = naiveGrobnerBasis(start);
          console.log("Iteration " + i + " complete");
            if (i == 9) {
            console.log(`Output basis size: ${basis.length}`);
            for (let i = 0; i < basis.length; i++) {
              console.log(`  [${i}]: ${basis[i].toString()}`);
            }
          }
        }
      }
      else {
        console.log("TypeScript generic finite coeff bitpacked exponent cyclic 5");
        const p1 = new Polynomial<IntModP, BitPackedExponents>([
          { coefficient: new IntModP(1), exponents: BitPackedExponents.fromArray([1, 0, 0, 0, 0, 0]) },
          { coefficient: new IntModP(1), exponents: BitPackedExponents.fromArray([0, 1, 0, 0, 0, 0]) },
          { coefficient: new IntModP(1), exponents: BitPackedExponents.fromArray([0, 0, 1, 0, 0, 0]) },
          { coefficient: new IntModP(1), exponents: BitPackedExponents.fromArray([0, 0, 0, 1, 0, 0]) },
          { coefficient: new IntModP(1), exponents: BitPackedExponents.fromArray([0, 0, 0, 0, 1, 0]) },
        ]);
        const p2 = new Polynomial<IntModP, BitPackedExponents>([
          { coefficient: new IntModP(1), exponents: BitPackedExponents.fromArray([1, 1, 0, 0, 0, 0]) },
          { coefficient: new IntModP(1), exponents: BitPackedExponents.fromArray([0, 1, 1, 0, 0, 0]) },
          { coefficient: new IntModP(1), exponents: BitPackedExponents.fromArray([0, 0, 1, 1, 0, 0]) },
          { coefficient: new IntModP(1), exponents: BitPackedExponents.fromArray([0, 0, 0, 1, 1, 0]) },
          { coefficient: new IntModP(1), exponents: BitPackedExponents.fromArray([1, 0, 0, 0, 1, 0]) },
        ]);
        const p3 = new Polynomial<IntModP, BitPackedExponents>([
          { coefficient: new IntModP(1), exponents: BitPackedExponents.fromArray([1, 1, 1, 0, 0, 0]) },
          { coefficient: new IntModP(1), exponents: BitPackedExponents.fromArray([0, 1, 1, 1, 0, 0]) },
          { coefficient: new IntModP(1), exponents: BitPackedExponents.fromArray([0, 0, 1, 1, 1, 0]) },
          { coefficient: new IntModP(1), exponents: BitPackedExponents.fromArray([1, 0, 0, 1, 1, 0]) },
          { coefficient: new IntModP(1), exponents: BitPackedExponents.fromArray([1, 1, 0, 0, 1, 0]) },
        ]);
        const p4 = new Polynomial<IntModP, BitPackedExponents>([
          { coefficient: new IntModP(1), exponents: BitPackedExponents.fromArray([1, 1, 1, 1, 0, 0]) },
          { coefficient: new IntModP(1), exponents: BitPackedExponents.fromArray([0, 1, 1, 1, 1, 0]) },
          { coefficient: new IntModP(1), exponents: BitPackedExponents.fromArray([1, 0, 1, 1, 1, 0]) },
          { coefficient: new IntModP(1), exponents: BitPackedExponents.fromArray([1, 1, 0, 1, 1, 0]) },
          { coefficient: new IntModP(1), exponents: BitPackedExponents.fromArray([1, 1, 1, 0, 1, 0]) },
        ]);
        const p5 = new Polynomial<IntModP, BitPackedExponents>([
          { coefficient: new IntModP(1), exponents: BitPackedExponents.fromArray([1, 1, 1, 1, 1, 0]) },
          { coefficient: new IntModP(modulus - 1), exponents: BitPackedExponents.fromArray([0, 0, 0, 0, 0, 0]) },
        ]);
        const start = [p1, p2, p3, p4, p5];
        for (let i = 0; i < 10; i++) {
          const basis = naiveGrobnerBasis(start);
          console.log("Iteration " + i + " complete");
            if (i == 9) {
            console.log(`Output basis size: ${basis.length}`);
            for (let i = 0; i < basis.length; i++) {
              console.log(`  [${i}]: ${basis[i].toString()}`);
            }
          }
        }
      }
    }
    else if (n == 6){
      if (vecType == 0) {
        console.log("TypeScript generic finite coeff vec exponent cyclic 6");
        const p1 = new Polynomial<IntModP, VecExponents>([
          { coefficient: new IntModP(1), exponents: new VecExponents([1, 0, 0, 0, 0, 0]) },
          { coefficient: new IntModP(1), exponents: new VecExponents([0, 1, 0, 0, 0, 0]) },
          { coefficient: new IntModP(1), exponents: new VecExponents([0, 0, 1, 0, 0, 0]) },
          { coefficient: new IntModP(1), exponents: new VecExponents([0, 0, 0, 1, 0, 0]) },
          { coefficient: new IntModP(1), exponents: new VecExponents([0, 0, 0, 0, 1, 0]) },
          { coefficient: new IntModP(1), exponents: new VecExponents([0, 0, 0, 0, 0, 1]) },
        ]);
        const p2 = new Polynomial<IntModP, VecExponents>([
          { coefficient: new IntModP(1), exponents: new VecExponents([1, 1, 0, 0, 0, 0]) },
          { coefficient: new IntModP(1), exponents: new VecExponents([0, 1, 1, 0, 0, 0]) },
          { coefficient: new IntModP(1), exponents: new VecExponents([0, 0, 1, 1, 0, 0]) },
          { coefficient: new IntModP(1), exponents: new VecExponents([0, 0, 0, 1, 1, 0]) },
          { coefficient: new IntModP(1), exponents: new VecExponents([0, 0, 0, 0, 1, 1]) },
          { coefficient: new IntModP(1), exponents: new VecExponents([1, 0, 0, 0, 0, 1]) },
        ]);
        const p3 = new Polynomial<IntModP, VecExponents>([
          { coefficient: new IntModP(1), exponents: new VecExponents([1, 1, 1, 0, 0, 0]) },
          { coefficient: new IntModP(1), exponents: new VecExponents([0, 1, 1, 1, 0, 0]) },
          { coefficient: new IntModP(1), exponents: new VecExponents([0, 0, 1, 1, 1, 0]) },
          { coefficient: new IntModP(1), exponents: new VecExponents([0, 0, 0, 1, 1, 1]) },
          { coefficient: new IntModP(1), exponents: new VecExponents([1, 0, 0, 0, 1, 1]) },
          { coefficient: new IntModP(1), exponents: new VecExponents([1, 1, 0, 0, 0, 1]) },
        ]);
        const p4 = new Polynomial<IntModP, VecExponents>([
          { coefficient: new IntModP(1), exponents: new VecExponents([1, 1, 1, 1, 0, 0]) },
          { coefficient: new IntModP(1), exponents: new VecExponents([0, 1, 1, 1, 1, 0]) },
          { coefficient: new IntModP(1), exponents: new VecExponents([0, 0, 1, 1, 1, 1]) },
          { coefficient: new IntModP(1), exponents: new VecExponents([1, 0, 0, 1, 1, 1]) },
          { coefficient: new IntModP(1), exponents: new VecExponents([1, 1, 0, 0, 1, 1]) },
          { coefficient: new IntModP(1), exponents: new VecExponents([1, 1, 1, 0, 0, 1]) },
        ]);
        const p5 = new Polynomial<IntModP, VecExponents>([
          { coefficient: new IntModP(1), exponents: new VecExponents([1, 1, 1, 1, 1, 0]) },
          { coefficient: new IntModP(1), exponents: new VecExponents([0, 1, 1, 1, 1, 1]) },
          { coefficient: new IntModP(1), exponents: new VecExponents([1, 0, 1, 1, 1, 1]) },
          { coefficient: new IntModP(1), exponents: new VecExponents([1, 1, 0, 1, 1, 1]) },
          { coefficient: new IntModP(1), exponents: new VecExponents([1, 1, 1, 0, 1, 1]) },
          { coefficient: new IntModP(1), exponents: new VecExponents([1, 1, 1, 1, 0, 1]) },
        ]);
        const p6 = new Polynomial<IntModP, VecExponents>([
          { coefficient: new IntModP(1), exponents: new VecExponents([1, 1, 1, 1, 1, 1]) },
          { coefficient: new IntModP(modulus - 1), exponents: new VecExponents([0, 0, 0, 0, 0, 0]) },
        ]);
        const start = [p1, p2, p3, p4, p5, p6];
        for (let i = 0; i < 10; i++) {
          const basis = naiveGrobnerBasis(start);
          console.log("Iteration " + i + " complete");
            if (i == 9) {
            console.log(`Output basis size: ${basis.length}`);
            for (let i = 0; i < basis.length; i++) {
              console.log(`  [${i}]: ${basis[i].toString()}`);
            }
          }
        }
      }
      else {
        console.log("TypeScript generic finite coeff bitpacked exponent cyclic 6");
        const p1 = new Polynomial<IntModP, BitPackedExponents>([
          { coefficient: new IntModP(1), exponents: BitPackedExponents.fromArray([1, 0, 0, 0, 0, 0]) },
          { coefficient: new IntModP(1), exponents: BitPackedExponents.fromArray([0, 1, 0, 0, 0, 0]) },
          { coefficient: new IntModP(1), exponents: BitPackedExponents.fromArray([0, 0, 1, 0, 0, 0]) },
          { coefficient: new IntModP(1), exponents: BitPackedExponents.fromArray([0, 0, 0, 1, 0, 0]) },
          { coefficient: new IntModP(1), exponents: BitPackedExponents.fromArray([0, 0, 0, 0, 1, 0]) },
          { coefficient: new IntModP(1), exponents: BitPackedExponents.fromArray([0, 0, 0, 0, 0, 1]) },
        ]);
        const p2 = new Polynomial<IntModP, BitPackedExponents>([
          { coefficient: new IntModP(1), exponents: BitPackedExponents.fromArray([1, 1, 0, 0, 0, 0]) },
          { coefficient: new IntModP(1), exponents: BitPackedExponents.fromArray([0, 1, 1, 0, 0, 0]) },
          { coefficient: new IntModP(1), exponents: BitPackedExponents.fromArray([0, 0, 1, 1, 0, 0]) },
          { coefficient: new IntModP(1), exponents: BitPackedExponents.fromArray([0, 0, 0, 1, 1, 0]) },
          { coefficient: new IntModP(1), exponents: BitPackedExponents.fromArray([0, 0, 0, 0, 1, 1]) },
          { coefficient: new IntModP(1), exponents: BitPackedExponents.fromArray([1, 0, 0, 0, 0, 1]) },
        ]);
        const p3 = new Polynomial<IntModP, BitPackedExponents>([
          { coefficient: new IntModP(1), exponents: BitPackedExponents.fromArray([1, 1, 1, 0, 0, 0]) },
          { coefficient: new IntModP(1), exponents: BitPackedExponents.fromArray([0, 1, 1, 1, 0, 0]) },
          { coefficient: new IntModP(1), exponents: BitPackedExponents.fromArray([0, 0, 1, 1, 1, 0]) },
          { coefficient: new IntModP(1), exponents: BitPackedExponents.fromArray([0, 0, 0, 1, 1, 1]) },
          { coefficient: new IntModP(1), exponents: BitPackedExponents.fromArray([1, 0, 0, 0, 1, 1]) },
          { coefficient: new IntModP(1), exponents: BitPackedExponents.fromArray([1, 1, 0, 0, 0, 1]) },
        ]);
        const p4 = new Polynomial<IntModP, BitPackedExponents>([
          { coefficient: new IntModP(1), exponents: BitPackedExponents.fromArray([1, 1, 1, 1, 0, 0]) },
          { coefficient: new IntModP(1), exponents: BitPackedExponents.fromArray([0, 1, 1, 1, 1, 0]) },
          { coefficient: new IntModP(1), exponents: BitPackedExponents.fromArray([0, 0, 1, 1, 1, 1]) },
          { coefficient: new IntModP(1), exponents: BitPackedExponents.fromArray([1, 0, 0, 1, 1, 1]) },
          { coefficient: new IntModP(1), exponents: BitPackedExponents.fromArray([1, 1, 0, 0, 1, 1]) },
          { coefficient: new IntModP(1), exponents: BitPackedExponents.fromArray([1, 1, 1, 0, 0, 1]) },
        ]);
        const p5 = new Polynomial<IntModP, BitPackedExponents>([
          { coefficient: new IntModP(1), exponents: BitPackedExponents.fromArray([1, 1, 1, 1, 1, 0]) },
          { coefficient: new IntModP(1), exponents: BitPackedExponents.fromArray([0, 1, 1, 1, 1, 1]) },
          { coefficient: new IntModP(1), exponents: BitPackedExponents.fromArray([1, 0, 1, 1, 1, 1]) },
          { coefficient: new IntModP(1), exponents: BitPackedExponents.fromArray([1, 1, 0, 1, 1, 1]) },
          { coefficient: new IntModP(1), exponents: BitPackedExponents.fromArray([1, 1, 1, 0, 1, 1]) },
          { coefficient: new IntModP(1), exponents: BitPackedExponents.fromArray([1, 1, 1, 1, 0, 1]) },
        ]);
        const p6 = new Polynomial<IntModP, BitPackedExponents>([
          { coefficient: new IntModP(1), exponents: BitPackedExponents.fromArray([1, 1, 1, 1, 1, 1]) },
          { coefficient: new IntModP(modulus - 1), exponents: BitPackedExponents.fromArray([0, 0, 0, 0, 0, 0]) },
        ]);
        const start = [p1, p2, p3, p4, p5, p6];
        for (let i = 0; i < 10; i++) {
          const basis = naiveGrobnerBasis(start);
          console.log("Iteration " + i + " complete");
            if (i == 9) {
            console.log(`Output basis size: ${basis.length}`);
            for (let i = 0; i < basis.length; i++) {
              console.log(`  [${i}]: ${basis[i].toString()}`);
            }
          }
        }
      }
    }
  }
}

main();
