// Specialized Grobner basis with bitpacked exponents (Rust-style)
// Supports up to 6 variables, 8 bits per variable, degree in top 16 bits

import { LCG } from "../helpers/lcg";

let modulus: number = 13;

export enum TermOrder {
  Lex,
  GrLex,
  RevLex,
}

let termOrder = TermOrder.Lex;

export type Term = {
  coefficient: number;
  exponents: bigint; // 64 bits: [63..48]=degree, [47..40]=e0, ... [7..0]=e5
}

function fromExponents(coefficient: number, exponents: number[]): Term {
    let packed = 0n;
    let degree = 0;
    for (let i = 0; i < 6; i++) {
      const exp = BigInt(exponents[i] ?? 0);
      packed |= exp << BigInt(40 - 8 * i);
      degree += Number(exp);
    }
    packed |= BigInt(degree) << 48n;
    return { coefficient, exponents: packed };
}

function unpack(term: Term): number[] {
  const exps: number[] = [];
  for (let i = 0; i < 6; i++) {
    exps.push(Number((term.exponents >> BigInt(40 - 8 * i)) & 0xFFn));
  }
  return exps;
}

function degree(term: Term): number {
  return Number((term.exponents >> 48n) & 0xFFFFn);
}

function compare(term: Term, other: Term): number {
  switch (termOrder) {
    case TermOrder.Lex: {
      const a = term.exponents & 0x0000FFFFFFFFFFFFn;
      const b = other.exponents & 0x0000FFFFFFFFFFFFn;
      return a > b ? 1 : a < b ? -1 : 0;
    }
    case TermOrder.GrLex: {
      const da = degree(term);
      const db = degree(other);
      if (da !== db) return da - db;
      const a = term.exponents & 0x0000FFFFFFFFFFFFn;
      const b = other.exponents & 0x0000FFFFFFFFFFFFn;
      return a > b ? 1 : a < b ? -1 : 0;
    }
    case TermOrder.RevLex: {
      const da = degree(term);
      const db = degree(other);
      if (da !== db) return da - db;
      const a = term.exponents & 0x0000FFFFFFFFFFFFn;
      const b = other.exponents & 0x0000FFFFFFFFFFFFn;
      return b > a ? 1 : b < a ? -1 : 0;
    }
    default:
      return 0;
  }
}

function canReduce(term: Term, divisor: Term): boolean {
  for (let i = 0; i < 6; i++) {
    const selfExp = Number((term.exponents >> BigInt(40 - 8 * i)) & 0xFFn);
    const divExp = Number((divisor.exponents >> BigInt(40 - 8 * i)) & 0xFFn);
    if (divExp > selfExp) return false;
  }
  return true;
}

function lcm(term: Term, other: Term): bigint {
  let lcmPacked = 0n;
  let degree = 0;
  for (let i = 0; i < 6; i++) {
    const a = Number((term.exponents >> BigInt(40 - 8 * i)) & 0xFFn);
    const b = Number((other.exponents >> BigInt(40 - 8 * i)) & 0xFFn);
    const lcmExp = Math.max(a, b);
    lcmPacked |= BigInt(lcmExp) << BigInt(40 - 8 * i);
    degree += lcmExp;
  }
  lcmPacked |= BigInt(degree) << 48n;
  return lcmPacked;
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

function polyKey(terms: Term[]): string {
  let out = "";
  for (let i = 0; i < terms.length; i++) {
    if (i > 0) out += ";";
    out += terms[i].coefficient + "|";
    out += terms[i].exponents + ", ";
  }
  return out;
}

export class Polynomial {
  terms: Term[];
  constructor(terms: Term[]) {
    this.terms = terms.filter(t => t.coefficient % modulus !== 0);
    this.terms = this.terms.map(t => ({
      coefficient: ((t.coefficient % modulus) + modulus) % modulus,
      exponents: t.exponents,
    }));
    this.terms.sort((a, b) => -compare(a,b));
  }

  add(other: Polynomial): Polynomial {
    const result: Term[] = this.terms.map(t => ({ coefficient: t.coefficient, exponents: t.exponents }));
    for (const t of other.terms) {
      let found = false;
      for (const rt of result) {
        if (rt.exponents === t.exponents) {
          rt.coefficient = (rt.coefficient + t.coefficient) % modulus;
          found = true;
          break;
        }
      }
      if (!found) result.push({ coefficient: t.coefficient, exponents: t.exponents });
    }
    return new Polynomial(result);
  }

  subtract(other: Polynomial): Polynomial {
    const result: Term[] = this.terms.map(t => ({ coefficient: t.coefficient, exponents: t.exponents }));
    for (const t of other.terms) {
      let found = false;
      for (const rt of result) {
        if (rt.exponents === t.exponents) {
          rt.coefficient = (rt.coefficient - t.coefficient + modulus) % modulus;
          found = true;
          break;
        }
      }
      if (!found) result.push({ coefficient: ((-t.coefficient + modulus) % modulus), exponents: t.exponents });
    }
    return new Polynomial(result);
  }

  multiplyByTerm(term: Term): Polynomial {
    const terms: Term[] = this.terms.map(t => ({
      coefficient: (t.coefficient * term.coefficient) % modulus,
      exponents: t.exponents + term.exponents,
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
        if (canReduce(lead, divLead)) {
          const coeff = (lead.coefficient * modInverse(divLead.coefficient, modulus)) % modulus;
          const exps = lead.exponents - divLead.exponents;
          const reductionTerm: Term = { coefficient: coeff, exponents: exps };
          const scaledDivisor = divisor.multiplyByTerm(reductionTerm);
          result = result.subtract(scaledDivisor);
          reduced = true;
          break;
        }
      }
      if (!reduced){ 
        if (result.terms.length === 0) break;
        remainder.push(result.terms[0]);
        result = new Polynomial(result.terms.slice(1));
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
      exponents: t.exponents,
    }));
    return new Polynomial(newTerms);
  }

  static sPolynomial(p1: Polynomial, p2: Polynomial): Polynomial {
    const lead1 = p1.terms[0];
    const lead2 = p2.terms[0];
    const lcmExps = lcm(lead1, lead2);
    const scale1: Term = { coefficient: lead2.coefficient, exponents: lcmExps - lead1.exponents };
    const scale2: Term = { coefficient: lead1.coefficient, exponents: lcmExps - lead2.exponents };
    const scaled1 = p1.multiplyByTerm(scale1);
    const scaled2 = p2.multiplyByTerm(scale2);
    return scaled1.subtract(scaled2);
  }
}

function generateCyclicPolynomials(n: number): Polynomial[] {
  // Generate cyclic-n polynomial system
  // For i = 0 to n-1: sum of all products of (i+1) consecutive variables = 0
  // Last equation: product of all n variables - 1 = 0
  
  const polys: Polynomial[] = [];
  
  for (let i = 0; i < n; i++) {
    const terms: Term[] = [];
    
    if (i === n - 1) {
      // Last equation: x0*x1*...*x(n-1) - 1
      const exps = Array(6).fill(0);
      for (let j = 0; j < n; j++) {
        exps[j] = 1;
      }
      terms.push(fromExponents(1, exps));
      terms.push(fromExponents((modulus - 1) % modulus, [0, 0, 0, 0, 0, 0])); // -1 mod m
    } else {
      // Equation i: sum of all products of (i+1) consecutive variables
      const productSize = i + 1;
      for (let start = 0; start < n; start++) {
        const exps = Array(6).fill(0);
        for (let j = 0; j < productSize; j++) {
          exps[(start + j) % n] = 1;
        }
        terms.push(fromExponents(1, exps));
      }
    }
    
    polys.push(new Polynomial(terms));
  }
  
  return polys;
}

export function naiveGrobnerBasis(polys: Polynomial[]): Polynomial[] {
  let basis = polys.map(p => new Polynomial(p.terms));
  const basisSet = new Set<string>();
   // Initialize basis set
    for (const poly of basis) {
      basisSet.add(polyKey(poly.terms));
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
      const key = polyKey(reduced.terms);
  
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
  const reducedSet = new Set<string>();
  for (const poly of basis) {
    const basisExcl = basis.filter(p => p !== poly);
    const reduced = poly.reduce(basisExcl);
    const key = polyKey(reduced.terms);
    if (reduced.terms.length > 0 && !reducedSet.has(key)) {
      reducedSet.add(key);
      reducedBasis.push(reduced.makeMonic());
    } 
  }
  return reducedBasis;
}

function main() {
  modulus = 7;
  const args = process.argv.slice(2);
  const n = args[0] ? parseInt(args[0], 10) : 4;
  // Test cyclic-4, cyclic-5, cyclic-6
  
  
  const polys = generateCyclicPolynomials(n);
  console.log(`Typescript specialized finite coeff bitpacked exponent cyclic ${n}`);
  
  for (let i = 0; i < 10; i++) {
    const basis = naiveGrobnerBasis(polys);
    console.log("Iteration", i, " complete");
      if (i === 9) {
      console.log(`Output basis size: ${basis.length}`);
      
      // Print final basis
      for (let i = 0; i < basis.length; i++) {
        const expsStr = basis[i].terms.map(t => {
          const exps = unpack(t);
          return `${t.coefficient}*${exps.slice(0, n).join('|')}`;
        }).join(' + ');
        console.log(`  [${i}]: ${expsStr}`);
      }
    }
  }

}

main();