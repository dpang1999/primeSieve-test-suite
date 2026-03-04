
// --- Finite field Grobner basis, grobner.ts style ---

import { LCG } from "../helpers/lcg";

let modulus: number = 13; // Global modulus for all terms, can be set from main

export type Term = {
  coefficient: number;
  exponents: number[];
};

export enum TermOrder {
  Lex,
  GrLex,
  RevLex,
}
let termOrder = TermOrder.Lex;

function compareExponents(a: number[], b: number[]): number {
  switch (termOrder) {
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
  constructor(terms: Term[]) {
    // Remove zero coefficients, sort by term order
    this.terms = terms.filter(t => t.coefficient % modulus !== 0);
    this.terms = this.terms.map(t => ({
      coefficient: ((t.coefficient % modulus) + modulus) % modulus,
      modulus: modulus,
      exponents: t.exponents.slice(),
    }));
    this.terms.sort((a, b) => -compareExponents(a.exponents, b.exponents));
  }

  add(other: Polynomial): Polynomial {
    const result = this.terms.map(t => ({ ...t }));
    for (const t of other.terms) {
      let found = false;
      for (const rt of result) {
        if (rt.exponents.length === t.exponents.length && rt.exponents.every((v, i) => v === t.exponents[i])) {
          rt.coefficient = (rt.coefficient + t.coefficient) % modulus;
          found = true;
          break;
        }
      }
      if (!found) result.push({ ...t });
    }
    return new Polynomial(result);
  }

  subtract(other: Polynomial): Polynomial {
    const result = this.terms.map(t => ({ ...t }));
    for (const t of other.terms) {
      let found = false;
      for (const rt of result) {
        if (rt.exponents.length === t.exponents.length && rt.exponents.every((v, i) => v === t.exponents[i])) {
          rt.coefficient = (rt.coefficient - t.coefficient + modulus) % modulus;
          found = true;
          break;
        }
      }
      if (!found) result.push({ coefficient: (-t.coefficient + modulus) % modulus, exponents: t.exponents.slice() });
    }
    return new Polynomial(result);
  }



  multiplyByTerm(term: Term): Polynomial {
    const terms = this.terms.map(t => ({
      coefficient: (t.coefficient * term.coefficient) % modulus,
      exponents: t.exponents.map((v, i) => v + term.exponents[i]),
    }));
    return new Polynomial(terms);
  }

  reduce(divisors: Polynomial[]): Polynomial {
    let result = new Polynomial(this.terms);
    const remainder: Term[] = [];

    while (true) {
      let reduced = false;
      for (const divisor of divisors) {
        if (result.terms.length === 0 || divisor.terms.length === 0) continue;
        const lead = result.terms[0];
        const divLead = divisor.terms[0];
        if (lead.exponents.length !== divLead.exponents.length) continue;
        if (lead.exponents.every((v, i) => v >= divLead.exponents[i])) {
          const coeff = (lead.coefficient * modInverse(divLead.coefficient, modulus)) % modulus;
          const exps = lead.exponents.map((v, i) => v - divLead.exponents[i]);
          const reductionTerm: Term = { coefficient: coeff, exponents: exps };
          const scaledDivisor = divisor.multiplyByTerm(reductionTerm);
          result = result.subtract(scaledDivisor);
          reduced = true;
          break;
        }
      }
      if (!reduced) 
      {
        if (result.terms.length === 0) break;
        remainder.push(result.terms[0]);
        result.terms = result.terms.slice(1);
      }
    }
    result.terms.push(...remainder);
    return new Polynomial(result.terms);
  }

  makeMonic(): Polynomial {
    if (this.terms.length === 0) return this;
    const leadCoeff = this.terms[0].coefficient;
    if (leadCoeff === 0) return this;
    const inv = modInverse(leadCoeff, modulus);
    const newTerms = this.terms.map(t => ({
      coefficient: (t.coefficient * inv) % modulus,
      exponents: t.exponents.slice(),
    }));
    return new Polynomial(newTerms);
  }

  static sPolynomial(p1: Polynomial, p2: Polynomial): Polynomial {
    const lead1 = p1.terms[0];
    const lead2 = p2.terms[0];
    const lcmExps = lead1.exponents.map((v, i) => Math.max(v, lead2.exponents[i]));
    const scale1 = lcmExps.map((v, i) => v - lead1.exponents[i]);
    const scale2 = lcmExps.map((v, i) => v - lead2.exponents[i]);
    const scaled1 = p1.multiplyByTerm({ coefficient: 1, exponents: scale1 });
    const scaled2 = p2.multiplyByTerm({ coefficient: 1, exponents: scale2 });
    return scaled1.subtract(scaled2);
  }
}

export function naiveGrobnerBasis(polys: Polynomial[]): Polynomial[] {
  let basis = polys.map(p => new Polynomial(p.terms));
  const basisSet = new Set<string>();
  
  // Initialize basis set
  for (const poly of basis) {
    basisSet.add(JSON.stringify(poly.terms));
  }

  // Initialize pairs: all (i, j) where i < j
  const pairs: Array<[number, number]> = [];
  for (let i = 0; i < basis.length; i++) {
    for (let j = i + 1; j < basis.length; j++) {
      pairs.push([i, j]);
    }
  }

  // Process pairs until none remain
  while (pairs.length > 0) {
    const [i, j] = pairs.shift()!;
    const sPoly = Polynomial.sPolynomial(basis[i], basis[j]);
    const reduced = sPoly.reduce(basis);
    const key = JSON.stringify(reduced.terms);

    if (reduced.terms.length > 0 && !basisSet.has(key)) {
      basisSet.add(key);
      const newIdx = basis.length;
      basis.push(reduced);
      
      // Add new pairs with the new polynomial
      for (let k = 0; k < newIdx; k++) {
        pairs.push([k, newIdx]);
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
      reducedBasis.push(reduced.makeMonic());
    }
  }
  return reducedBasis;
}


function main() {
  // let mode == 0 for testing
  const args = process.argv.slice(2);
  const mode = 0;
  if (mode != 0) {
    // arg1 = number of polynomials
    // arg2 = term order (0=Lex, 1=GrLex, 2=RevLex)
    // arg3 = modulus
    
    const numPolys = parseInt(args[0] || '3', 10);
    const orderArg = parseInt(args[1] || '0', 10);
    modulus = parseInt(args[2] || '13', 10);
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
        terms.push({ coefficient, exponents });
        }
        polys.push(new Polynomial(terms));
    }
    const basis = naiveGrobnerBasis(polys);
    console.log('Computed Grobner Basis:');
    for (const poly of basis) {
        for (const term of poly.terms) {
        process.stdout.write(`${term.coefficient}*x^${term.exponents[0]}y^${term.exponents[1]}z^${term.exponents[2]} + `);
        }
        console.log();
    }
  }
  else {
    const n = args[0] ? parseInt(args[0], 10) : 4;
    if ([4, 5, 6, 7].includes(n)) {
      modulus = 7;
      const polys: Polynomial[] = [];
      // f1: x0 + x1 + ... + x_{n-1}
      const f1_terms: Term[] = [];
      for (let i = 0; i < n; i++) {
        const exponents = Array(n).fill(0);
        exponents[i] = 1;
        f1_terms.push({ coefficient: 1, exponents });
      }
      polys.push(new Polynomial(f1_terms));
      // f2: x0x1 + x1x2 + ... + x_{n-1}x0
      const f2_terms: Term[] = [];
      for (let i = 0; i < n; i++) {
        const exponents = Array(n).fill(0);
        exponents[i] = 1;
        exponents[(i + 1) % n] = 1;
        f2_terms.push({ coefficient: 1, exponents });
      }
      polys.push(new Polynomial(f2_terms));
      // f3: x0x1x2 + x1x2x3 + ... + x_{n-1}x0x1
      const f3_terms: Term[] = [];
      for (let i = 0; i < n; i++) {
        const exponents = Array(n).fill(0);
        exponents[i] = 1;
        exponents[(i + 1) % n] = 1;
        exponents[(i + 2) % n] = 1;
        f3_terms.push({ coefficient: 1, exponents });
      }
      polys.push(new Polynomial(f3_terms));
      // Continue for f4, f5, ..., f_{n-1}
      for (let k = 4; k < n; k++) {
        const fk_terms: Term[] = [];
        for (let i = 0; i < n; i++) {
          const exponents = Array(n).fill(0);
          for (let j = 0; j < k; j++) {
            exponents[(i + j) % n] = 1;
          }
          fk_terms.push({ coefficient: 1, exponents });
        }
        polys.push(new Polynomial(fk_terms));
      }
      // fn: x0x1...x_{n-1} - 1
      const fn_terms: Term[] = [];
      fn_terms.push({ coefficient: 1, exponents: Array(n).fill(1) });
      fn_terms.push({ coefficient: modulus - 1, exponents: Array(n).fill(0) }); // -1 mod modulus
      polys.push(new Polynomial(fn_terms));

      // print polys
      console.log(`Input Polynomials for Cyclic-${n}:`);
      for (const poly of polys) {
        for (const term of poly.terms) {
          process.stdout.write(`${term.coefficient}*`);
          for (let i = 0; i < term.exponents.length; i++) {
            process.stdout.write(`x${i}^${term.exponents[i]}`);
            if (i < term.exponents.length - 1) process.stdout.write(" ");
          }
          process.stdout.write(" + ");
        }
        console.log();
      }

      const basis = naiveGrobnerBasis(polys);
      console.log(`Final Grobner Basis for Cyclic-${n}:`);
      console.log(`Number of basis elements: ${basis.length}`);
      for (const poly of basis) {
        for (const term of poly.terms) {
          process.stdout.write(`${term.coefficient}*`);
          for (let i = 0; i < term.exponents.length; i++) {
            process.stdout.write(`x${i}^${term.exponents[i]}`);
            if (i < term.exponents.length - 1) process.stdout.write(" ");
          }
          process.stdout.write(" + ");
        }
        console.log();
      }
    }
  }
}

if (require.main === module) {
  main();
}
