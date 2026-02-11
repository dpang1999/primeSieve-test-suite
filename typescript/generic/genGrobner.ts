import { LCG } from "../helpers/lcg";
import { IField } from "./iField";
import { IExponent } from "./iExponent";
export class Term<C extends IField<C>, E extends IExponent<E>> {
  coefficient: C;
  exponents: E;
  
  constructor(coefficient: C, exponents: E) {
    this.coefficient = coefficient;
    this.exponents = exponents;
  }
  
  compare(other: Term<C, E>, order: TermOrder): number {
    switch (order) {
      case TermOrder.Lex:
        return this.exponents.lexCompare(other.exponents);
      case TermOrder.GrLex: {
        const degA = this.exponents.deg();
        const degB = other.exponents.deg();
        if (degA !== degB) return degA - degB;
        return this.exponents.lexCompare(other.exponents);
      }
      case TermOrder.RevLex: {
        const degA = this.exponents.deg();
        const degB = other.exponents.deg();
        if (degA !== degB) return degA - degB;
        else { return other.exponents.lexCompare(this.exponents); }
      }
      default:
        return 0;
    }
  }

  can_reduce(other: Term<C, E>): boolean {
    return this.exponents.canReduce(other.exponents);
  }

  lcm(other: Term<C, E>): E {
    return this.exponents.lcm(other.exponents);
  }

};

export enum TermOrder {
  Lex,
  GrLex,
  RevLex,
}



class Polynomial<C extends IField<C>, E extends IExponent<E>> {
  terms: Term<C, E>[];
  constructor(terms: Term<C, E>[], order: TermOrder) {
    // Remove near-zero coefficients, round, and sort
    this.terms = terms
      .filter(t => Math.abs(t.coefficient.coerce_to_number()) > 1e-2)
      .sort((a, b) => b.compare(a, order));
  }

  add(other: Polynomial<C, E>, order: TermOrder): Polynomial<C, E> {
    const results = this.terms.map(t => new Term<C, E>(t.coefficient.copy(), t.exponents));
    for (const t of other.terms) {
      let found = false;
      for (const rt of results) {
        if (rt.exponents.equals(t.exponents)) {
          rt.coefficient = rt.coefficient.a(t.coefficient);
          found = true;
          break;
        }
      }
      if (!found) this.terms.push(new Term<C, E>(t.coefficient, t.exponents));
    }
    return new Polynomial<C, E>(this.terms, order);
  }

  subtract(other: Polynomial<C, E>, order: TermOrder): Polynomial<C, E> {
    const results = this.terms.map(t => new Term<C, E>(t.coefficient.copy(), t.exponents));
    for (const t of other.terms) {
      let found = false;
      for (const rt of results) {
        if (rt.exponents.equals(t.exponents)) {
          rt.coefficient = rt.coefficient.s(t.coefficient);
          found = true;
          break;
        }
      }
      if (!found) results.push(new Term<C, E>(t.coefficient.zero().s(t.coefficient), t.exponents));
    }
    return new Polynomial(results, order);
  }

  multiplyByTerm(term: Term<C, E>, order: TermOrder): Polynomial<C, E> {
    const terms = this.terms.map(t => new Term<C, E>(
      t.coefficient.m(term.coefficient),
      t.exponents.add(term.exponents),
    ));
    return new Polynomial(terms, order);
  }

  reduce(divisors: Polynomial<C, E>[], order: TermOrder): Polynomial<C, E> {
    let result = new Polynomial(this.terms, order);
    while (true) {
      let reduced = false;
      for (const divisor of divisors) {
        if (result.terms.length === 0 || divisor.terms.length === 0) continue;
        const lead = result.terms[0];
        const divLead = divisor.terms[0];
        if (lead.can_reduce(divLead)) {
          const coeff = lead.coefficient.d(divLead.coefficient);
          const exps = lead.exponents.sub(divLead.exponents);
          const reductionTerm: Term<C,E> = new Term<C,E>(coeff, exps);
          const scaledDivisor = divisor.multiplyByTerm(reductionTerm, order);
          result = result.subtract(scaledDivisor, order);
          reduced = true;
          break;
        }
      }
      if (!reduced) break;
    }
    return new Polynomial(result.terms, order);
  }

  static sPolynomial<C extends IField<C>, E extends IExponent<E>>(p1: Polynomial<C, E>, p2: Polynomial<C, E>, order: TermOrder): Polynomial<C, E> {
    const lead1 = p1.terms[0];
    const lead2 = p2.terms[0];
    const lcmExps = lead1.lcm(lead2);
    const scale1 = lcmExps.sub(lead1.exponents);
    const scale2 = lcmExps.sub(lead2.exponents);
    const scaled1 = p1.multiplyByTerm(new Term<C,E>(lead1.coefficient.one(), scale1), order);
    const scaled2 = p2.multiplyByTerm(new Term<C,E>(lead2.coefficient.one(), scale2), order);
    return scaled1.subtract(scaled2, order);
  }

toString() {
    return this.terms.map(term => term.coefficient.toString() + "*" + term.exponents.toString()).join(" + ");
  }
}

export function naiveGrobnerBasis<C extends IField<C>, E extends IExponent<E>>(polys: Polynomial<C, E>[], order: TermOrder): Polynomial<C, E>[] {
  let basis = polys.map(p => new Polynomial(p.terms, order));
  const basisSet = new Set<string>();
  let added = true;
  while (added) {
    added = false;
    const n = basis.length;
    //console.log("Basis length:" + n);
    for (let i = 0; i < n; i++) {
      for (let j = i + 1; j < n; j++) {
        const sPoly = Polynomial.sPolynomial(basis[i], basis[j], order);
        const reduced = sPoly.reduce(basis, order);
        const key = JSON.stringify(reduced.terms);
        if (reduced.terms.length > 0 && !basisSet.has(key)) {
          //console.log({i}, {j});
          basisSet.add(key);
          basis.push(reduced);
          added = true;
        }
      }
    }
  }

  /*console.log("Basis before reduction:");
  for (const poly of basis) {
    let str = poly.terms.map(term => term.coefficient.toString() + "*" + term.exponents.toString()).join(" + ");
    console.log(str);
  }*/


  // Reduce basis by self
  const reducedBasis: Polynomial<C, E>[] = [];
  for (const poly of basis) {
    const basisExcl = basis.filter(p => p !== poly);
    const reduced = poly.reduce(basisExcl, order);
    const key = JSON.stringify(reduced.terms);
    if (reduced.terms.length > 0 && !reducedBasis.some(rb => JSON.stringify(rb.terms) === key)) {
      reducedBasis.push(reduced);
    }
  }
  return reducedBasis;
}




// --- New main for Grobner basis testing, modeled after Java implementation ---
import { DoubleField } from "./doubleField";
import { VecExponents } from "./vecExponents";
import { IntModP } from "./intModP";

function main() {
  // args: numPolys, numTerms, coeffType, expType, orderArg
  const args = process.argv.slice(2);
  const numPolys = parseInt(args[0] || '3', 10);
  const numTerms = parseInt(args[1] || '3', 10);
  const coeffType = parseInt(args[2] || '0', 10); // 0 = DoubleField
  const expType = parseInt(args[3] || '0', 10);   // 0 = VecExponent
  const orderArg = parseInt(args[4] || '0', 10);
  const modulus = parseInt(args[5] || '13', 10); // for IntModP if needed
  IntModP.setModulus(modulus);
  let order: TermOrder;
  switch (orderArg) {
    case 0: order = TermOrder.Lex; break;
    case 1: order = TermOrder.GrLex; break;
    case 2: order = TermOrder.RevLex; break;
    default: order = TermOrder.Lex;
  }
  console.log(`Using term order: ${TermOrder[order]}`);
  const rand = new LCG(12345, 1345, 16645, 1013904);

  if (coeffType === 0 && expType === 0) {
    // DoubleField + VecExponent
    const inputBasis: Polynomial<DoubleField, VecExponents>[] = [];
    for (let i = 0; i < numPolys; i++) {
      const terms: Term<DoubleField, VecExponents>[] = [];
      for (let j = 0; j < numTerms; j++) {
        const coefficient = new DoubleField(rand.nextDouble());
        const exponents = new VecExponents([
          rand.nextInt() % 4,
          rand.nextInt() % 4,
          rand.nextInt() % 4
        ]);
        terms.push(new Term(coefficient, exponents));
      }
      inputBasis.push(new Polynomial(terms, order));
    }
    console.log("Input Polynomials:");
    for (const poly of inputBasis) {
      let str = poly.terms.map(term => term.coefficient.toString() + "*" + term.exponents.toString()).join(" + ");
      console.log(str);
    }

    /*const temp = Polynomial.sPolynomial(inputBasis[1], inputBasis[2], order);
    console.log("S-Polynomial of first two input polynomials:");
    // print inputBasis[1] and inputBasis[2]
    console.log(inputBasis[1].terms.map(term => term.coefficient.toString() + "*" + term.exponents.toString()).join(" + "));
    console.log(inputBasis[2].terms.map(term => term.coefficient.toString() + "*" + term.exponents.toString()).join(" + "));
    let sPolyStr = temp.terms.map(term => term.coefficient.toString() + "*" + term.exponents.toString()).join(" + ");
    console.log(sPolyStr);*/

    const basis = naiveGrobnerBasis(inputBasis, order);
    console.log("Computed Grobner Basis Polynomials:");
    for (const poly of basis) {
      let str = poly.terms.map(term => term.coefficient.toString() + "*" + term.exponents.toString()).join(" + ");
      console.log(str);
    }


  } else {
    console.log("Only DoubleField + VecExponent supported in this demo main.");
  }
}

if (require.main === module) {
  main();
}
