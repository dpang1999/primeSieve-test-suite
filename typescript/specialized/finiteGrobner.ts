
// --- Finite field Grobner basis, grobner.ts style ---

import { LCG } from "../helpers/lcg";

export type Term = {
  coefficient: number;
  modulus: number;
  exponents: number[];
};

export enum TermOrder {
  Lex,
  GrLex,
  RevLex,
}

function compareExponents(a: number[], b: number[], order: TermOrder): number {
  switch (order) {
    case TermOrder.Lex:
      for (let i = 0; i < a.length; i++) {
        if (a[i] !== b[i]) return a[i] - b[i];
      }
      return 0;
    case TermOrder.GrLex: {
      const degA = a.reduce((s, x) => s + x, 0);
      const degB = b.reduce((s, x) => s + x, 0);
      if (degA !== degB) return degA - degB;
      for (let i = 0; i < a.length; i++) {
        if (a[i] !== b[i]) return a[i] - b[i];
      }
      return 0;
    }
    case TermOrder.RevLex: {
      const degA = a.reduce((s, x) => s + x, 0);
      const degB = b.reduce((s, x) => s + x, 0);
      if (degA !== degB) return degA - degB;
      for (let i = a.length - 1; i >= 0; i--) {
        if (a[i] !== b[i]) return a[i] - b[i];
      }
      return 0;
    }
    default:
      return 0;
  }
}

function modInverse(a: number, m: number): number {
  let m0 = m, x0 = 0, x1 = 1;
  if (m === 1) return 0;
  a = ((a % m) + m) % m;
  while (a > 1) {
    const q = Math.floor(a / m);
    [a, m] = [m, a % m];
    [x0, x1] = [x1 - q * x0, x0];
  }
  if (x1 < 0) x1 += m0;
  return x1;
}

export class Polynomial {
  terms: Term[];
  order: TermOrder;
  constructor(terms: Term[], order: TermOrder) {
    this.order = order;
    // Remove zero coefficients, sort by term order
    this.terms = terms.filter(t => t.coefficient % t.modulus !== 0);
    this.terms = this.terms.map(t => ({
      coefficient: ((t.coefficient % t.modulus) + t.modulus) % t.modulus,
      modulus: t.modulus,
      exponents: t.exponents.slice(),
    }));
    this.terms.sort((a, b) => -compareExponents(a.exponents, b.exponents, order));
  }

  add(other: Polynomial): Polynomial {
    const result = this.terms.map(t => ({ ...t }));
    for (const t of other.terms) {
      let found = false;
      for (const rt of result) {
        if (rt.exponents.length === t.exponents.length && rt.exponents.every((v, i) => v === t.exponents[i])) {
          rt.coefficient = (rt.coefficient + t.coefficient) % t.modulus;
          found = true;
          break;
        }
      }
      if (!found) result.push({ ...t });
    }
    return new Polynomial(result, this.order);
  }

  subtract(other: Polynomial): Polynomial {
    const result = this.terms.map(t => ({ ...t }));
    for (const t of other.terms) {
      let found = false;
      for (const rt of result) {
        if (rt.exponents.length === t.exponents.length && rt.exponents.every((v, i) => v === t.exponents[i])) {
          rt.coefficient = (rt.coefficient - t.coefficient + t.modulus) % t.modulus;
          found = true;
          break;
        }
      }
      if (!found) result.push({ coefficient: (-t.coefficient + t.modulus) % t.modulus, modulus: t.modulus, exponents: t.exponents.slice() });
    }
    return new Polynomial(result, this.order);
  }

  multiplyByTerm(term: Term): Polynomial {
    const terms = this.terms.map(t => ({
      coefficient: (t.coefficient * term.coefficient) % term.modulus,
      modulus: term.modulus,
      exponents: t.exponents.map((v, i) => v + term.exponents[i]),
    }));
    return new Polynomial(terms, this.order);
  }

  reduce(divisors: Polynomial[]): Polynomial {
    let result = new Polynomial(this.terms, this.order);
    while (true) {
      let reduced = false;
      for (const divisor of divisors) {
        if (result.terms.length === 0 || divisor.terms.length === 0) continue;
        const lead = result.terms[0];
        const divLead = divisor.terms[0];
        if (lead.exponents.length !== divLead.exponents.length) continue;
        if (lead.exponents.every((v, i) => v >= divLead.exponents[i])) {
          const coeff = (lead.coefficient * modInverse(divLead.coefficient, lead.modulus)) % lead.modulus;
          const exps = lead.exponents.map((v, i) => v - divLead.exponents[i]);
          const reductionTerm: Term = { coefficient: coeff, modulus: lead.modulus, exponents: exps };
          const scaledDivisor = divisor.multiplyByTerm(reductionTerm);
          result = result.subtract(scaledDivisor);
          reduced = true;
          break;
        }
      }
      if (!reduced) break;
    }
    return new Polynomial(result.terms, this.order);
  }

  static sPolynomial(p1: Polynomial, p2: Polynomial): Polynomial {
    const lead1 = p1.terms[0];
    const lead2 = p2.terms[0];
    const lcmExps = lead1.exponents.map((v, i) => Math.max(v, lead2.exponents[i]));
    const scale1 = lcmExps.map((v, i) => v - lead1.exponents[i]);
    const scale2 = lcmExps.map((v, i) => v - lead2.exponents[i]);
    const scaled1 = p1.multiplyByTerm({ coefficient: 1, modulus: lead1.modulus, exponents: scale1 });
    const scaled2 = p2.multiplyByTerm({ coefficient: 1, modulus: lead2.modulus, exponents: scale2 });
    return scaled1.subtract(scaled2);
  }
}

export function naiveGrobnerBasis(polys: Polynomial[], order: TermOrder): Polynomial[] {
  let basis = polys.map(p => new Polynomial(p.terms, order));
  const basisSet = new Set<string>();
  let added = true;
  while (added) {
    added = false;
    const n = basis.length;
    for (let i = 0; i < n; i++) {
      for (let j = i + 1; j < n; j++) {
        const sPoly = Polynomial.sPolynomial(basis[i], basis[j]);
        const reduced = sPoly.reduce(basis);
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
  const reducedBasis: Polynomial[] = [];
  for (const poly of basis) {
    const basisExcl = basis.filter(p => p !== poly);
    const reduced = poly.reduce(basisExcl);
    const key = JSON.stringify(reduced.terms);
    if (reduced.terms.length > 0 && !reducedBasis.some(rb => JSON.stringify(rb.terms) === key)) {
      reducedBasis.push(reduced);
    }
  }
  return reducedBasis;
}


function main() {
  // let mode == 0 for testing
  const mode = 0;
  if (mode != 0) {
    // arg1 = number of polynomials
    // arg2 = term order (0=Lex, 1=GrLex, 2=RevLex)
    // arg3 = modulus
    const args = process.argv.slice(2);
    const numPolys = parseInt(args[0] || '3', 10);
    const orderArg = parseInt(args[1] || '0', 10);
    const modulus = parseInt(args[2] || '13', 10);
    let order: TermOrder;
    switch (orderArg) {
      case 0:
        order = TermOrder.Lex;
        break;
      case 1:
        order = TermOrder.GrLex;
        break;
      case 2:
        order = TermOrder.RevLex;
        break;
      default:
        order = TermOrder.Lex;
    }
    const rand = new LCG(12345, 1345, 16645, 1013904);
    const polys: Polynomial[] = [];
    for (let i = 0; i < numPolys; ++i) {
        const terms: Term[] = [];
        for (let j = 0; j < 3; ++j) {
        const coefficient = rand.nextDouble() * 2 - 1;
        const exponents = [0, 0, 0].map(() => rand.nextInt() % 4);
        terms.push({ coefficient, modulus, exponents });
        }
        polys.push(new Polynomial(terms, order));
    }
    const basis = naiveGrobnerBasis(polys, order);
    console.log('Computed Grobner Basis:');
    for (const poly of basis) {
        for (const term of poly.terms) {
        process.stdout.write(`${term.coefficient}*x^${term.exponents[0]}y^${term.exponents[1]}z^${term.exponents[2]} + `);
        }
        console.log();
    }
  }
  else {
    // Example: x^3 + y^3 + z^3
    const modulus = 13;
    const p1 = new Polynomial([
      { coefficient: 1, modulus: 13, exponents: [3, 0, 0] },
      { coefficient: 1, modulus: 13, exponents: [0, 3, 0] },
      { coefficient: 1, modulus: 13, exponents: [0, 0, 3] },
    ], TermOrder.Lex);
    // xy + yz + xz
    const p2 = new Polynomial([
      { coefficient: 1, modulus: 13, exponents: [1, 1, 0] },
      { coefficient: 1, modulus: 13, exponents: [0, 1, 1] },
      { coefficient: 1, modulus: 13, exponents: [1, 0, 1] },
    ], TermOrder.Lex);
    // x+y+z
    const p = new Polynomial([
      { coefficient: 1, modulus: 13, exponents: [1, 0, 0] },
      { coefficient: 1, modulus: 13, exponents: [0, 1, 0] },
      { coefficient: 1, modulus: 13, exponents: [0, 0, 1] },
    ], TermOrder.Lex);
    const basis = naiveGrobnerBasis([p1, p2, p], TermOrder.Lex);
    console.log('Final Grobner Basis:');
    for (const poly of basis) {
      for (const term of poly.terms) {
        process.stdout.write(`${term.coefficient}*x^${term.exponents[0]}y^${term.exponents[1]}z^${term.exponents[2]} + `);
      }
      console.log();
    }
  }
}

if (require.main === module) {
  main();
}
