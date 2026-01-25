// Specialized Grobner basis with bitpacked exponents (Rust-style)
// Supports up to 6 variables, 8 bits per variable, degree in top 16 bits

import { LCG } from "../helpers/lcg";

export enum TermOrder {
  Lex,
  GrLex,
  RevLex,
}

export class Term {
  coefficient: number;
  exponents: bigint; // 64 bits: [63..48]=degree, [47..40]=e0, ... [7..0]=e5

  constructor(coefficient: number, exponents: bigint) {
    this.coefficient = coefficient;
    this.exponents = exponents;
  }

  static fromExponents(coefficient: number, exponents: number[]): Term {
    let packed = 0n;
    let degree = 0;
    for (let i = 0; i < 6; i++) {
      const exp = BigInt(exponents[i] ?? 0);
      packed |= exp << BigInt(40 - 8 * i);
      degree += Number(exp);
    }
    packed |= BigInt(degree) << 48n;
    return new Term(coefficient, packed);
  }

  unpack(): number[] {
    const exps: number[] = [];
    for (let i = 0; i < 6; i++) {
      exps.push(Number((this.exponents >> BigInt(40 - 8 * i)) & 0xFFn));
    }
    return exps;
  }

  degree(): number {
    return Number((this.exponents >> 48n) & 0xFFFFn);
  }

  compare(other: Term, order: TermOrder): number {
    switch (order) {
      case TermOrder.Lex: {
        const a = this.exponents & 0x0000FFFFFFFFFFFFn;
        const b = other.exponents & 0x0000FFFFFFFFFFFFn;
        return a > b ? 1 : a < b ? -1 : 0;
      }
      case TermOrder.GrLex: {
        const da = this.degree();
        const db = other.degree();
        if (da !== db) return da - db;
        const a = this.exponents & 0x0000FFFFFFFFFFFFn;
        const b = other.exponents & 0x0000FFFFFFFFFFFFn;
        return a > b ? 1 : a < b ? -1 : 0;
      }
      case TermOrder.RevLex: {
        const da = this.degree();
        const db = other.degree();
        if (da !== db) return da - db;
        const a = this.exponents & 0x0000FFFFFFFFFFFFn;
        const b = other.exponents & 0x0000FFFFFFFFFFFFn;
        return b > a ? 1 : b < a ? -1 : 0;
      }
      default:
        return 0;
    }
  }

  canReduce(divisor: Term): boolean {
    for (let i = 0; i < 6; i++) {
      const selfExp = Number((this.exponents >> BigInt(40 - 8 * i)) & 0xFFn);
      const divExp = Number((divisor.exponents >> BigInt(40 - 8 * i)) & 0xFFn);
      if (divExp > selfExp) return false;
    }
    return true;
  }

  lcm(other: Term): bigint {
    let lcmPacked = 0n;
    let degree = 0;
    for (let i = 0; i < 6; i++) {
      const a = Number((this.exponents >> BigInt(40 - 8 * i)) & 0xFFn);
      const b = Number((other.exponents >> BigInt(40 - 8 * i)) & 0xFFn);
      const lcmExp = Math.max(a, b);
      lcmPacked |= BigInt(lcmExp) << BigInt(40 - 8 * i);
      degree += lcmExp;
    }
    lcmPacked |= BigInt(degree) << 48n;
    return lcmPacked;
  }
/*
  addPacked(other: Term): bigint {
    let packed = 0n;
    let degree = 0;
    for (let i = 0; i < 6; i++) {
      const a = Number((this.exponents >> BigInt(40 - 8 * i)) & 0xFFn);
      const b = Number((other.exponents >> BigInt(40 - 8 * i)) & 0xFFn);
      const sum = a + b;
      packed |= BigInt(sum) << BigInt(40 - 8 * i);
      degree += sum;
    }
    packed |= BigInt(degree) << 48n;
    return packed;
  }

  subPacked(other: Term): bigint {
    let packed = 0n;
    let degree = 0;
    for (let i = 0; i < 6; i++) {
      const a = Number((this.exponents >> BigInt(40 - 8 * i)) & 0xFFn);
      const b = Number((other.exponents >> BigInt(40 - 8 * i)) & 0xFFn);
      const diff = a - b;
      packed |= BigInt(diff) << BigInt(40 - 8 * i);
      degree += diff;
    }
    packed |= BigInt(degree) << 48n;
    return packed;
  }*/
}

export class Polynomial {
  terms: Term[];
  constructor(terms: Term[], order: TermOrder) {
    this.terms = terms
      .filter(t => Math.abs(t.coefficient) > 1e-2)
      .map(t => new Term(Math.round(t.coefficient * 1e5) / 1e5, t.exponents));
    this.terms.sort((a, b) => -a.compare(b, order));
  }

  add(other: Polynomial, order: TermOrder): Polynomial {
    const result = this.terms.map(t => new Term(t.coefficient, t.exponents));
    for (const t of other.terms) {
      let found = false;
      for (const rt of result) {
        if (rt.exponents === t.exponents) {
          rt.coefficient += t.coefficient;
          found = true;
          break;
        }
      }
      if (!found) result.push(new Term(t.coefficient, t.exponents));
    }
    return new Polynomial(result, order);
  }

  subtract(other: Polynomial, order: TermOrder): Polynomial {
    const result = this.terms.map(t => new Term(t.coefficient, t.exponents));
    for (const t of other.terms) {
      let found = false;
      for (const rt of result) {
        if (rt.exponents === t.exponents) {
          rt.coefficient -= t.coefficient;
          found = true;
          break;
        }
      }
      if (!found) result.push(new Term(-t.coefficient, t.exponents));
    }
    return new Polynomial(result, order);
  }

  multiplyByTerm(term: Term, order: TermOrder): Polynomial {
    const terms = this.terms.map(t => new Term(
      t.coefficient * term.coefficient,
      t.exponents + term.exponents
    ));
    return new Polynomial(terms, order);
  }

  reduce(divisors: Polynomial[], order: TermOrder): Polynomial {
    let result = new Polynomial(this.terms, order);
    while (true) {
      let reduced = false;
      for (const divisor of divisors) {
        if (result.terms.length === 0 || divisor.terms.length === 0) continue;
        const lead = result.terms[0];
        const divLead = divisor.terms[0];
        if (lead.canReduce(divLead)) {
          const coeff = lead.coefficient / divLead.coefficient;
          const exps = lead.exponents - divLead.exponents;
          const reductionTerm = new Term(coeff, exps);
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

  static sPolynomial(p1: Polynomial, p2: Polynomial, order: TermOrder): Polynomial {
    const lead1 = p1.terms[0];
    const lead2 = p2.terms[0];
    const lcmExps = lead1.lcm(lead2);
    const scale1 = new Term(1, lcmExps - lead1.exponents);
    const scale2 = new Term(1, lcmExps - lead2.exponents);
    const scaled1 = p1.multiplyByTerm(scale1, order);
    const scaled2 = p2.multiplyByTerm(scale2, order);
    return scaled1.subtract(scaled2, order);
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
        const sPoly = Polynomial.sPolynomial(basis[i], basis[j], order);
        const reduced = sPoly.reduce(basis, order);
        const key = JSON.stringify(reduced.terms.map(t => [t.coefficient, t.exponents.toString()]));
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
    const reduced = poly.reduce(basisExcl, order);
    const key = JSON.stringify(reduced.terms.map(t => [t.coefficient, t.exponents.toString()]));
    if (reduced.terms.length > 0 && !reducedBasis.some(rb => JSON.stringify(rb.terms.map(t => [t.coefficient, t.exponents.toString()])) === key)) {
      reducedBasis.push(reduced);
    }
  }
  return reducedBasis;
}

function main() {
  // let mode = 0 be for testing
  const mode = 0;
  if (mode !== 0) {
    // arg1 = number of polynomials
    // arg2 = term order (0=Lex, 1=GrLex, 2=RevLex)
    const args = process.argv.slice(2);
    const numPolys = parseInt(args[0] || '3', 10);
    const orderArg = parseInt(args[1] || '0', 10);
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
    for (let i = 0; i < numPolys; i++) {
      const terms: Term[] = [];
      for (let j = 0; j < 3; j++) {
        const coeff =   rand.nextDouble() * 2 - 1;
        const exps = Array.from({ length: 6 }, () => (rand.nextInt() % 4));
        terms.push(Term.fromExponents(coeff, exps));
      }
      polys.push(new Polynomial(terms, order));
    }
    const basis = naiveGrobnerBasis(polys, order);
    console.log(basis.length);
    //console.log('Computed Grobner Basis:', basis);
    return;
  } else {
    const order = TermOrder.Lex;
    // Example: x^2*y + y^2*z + z^2*x, x*y*z - 1, x+y+z
    const p1 = new Polynomial([
      Term.fromExponents(1, [2, 1, 0, 0, 0, 0]),
      Term.fromExponents(1, [0, 2, 1, 0, 0, 0]),
      Term.fromExponents(1, [1, 0, 2, 0, 0, 0]),
    ], order);
    const p2 = new Polynomial([
      Term.fromExponents(1, [1, 1, 1, 0, 0, 0]),
      Term.fromExponents(-1, [0, 0, 0, 0, 0, 0]),
    ], order);
    const p3 = new Polynomial([
      Term.fromExponents(1, [1, 0, 0, 0, 0, 0]),
      Term.fromExponents(1, [0, 1, 0, 0, 0, 0]),
      Term.fromExponents(1, [0, 0, 1, 0, 0, 0]),
    ], order);
    const basis = naiveGrobnerBasis([p1, p2, p3], order);
    console.log('Final Grobner Basis:');
    for (const poly of basis) {
      console.log(poly.terms.map(t => ({ coeff: t.coefficient, exps: t.unpack() })));
    }
  }
}

if (require.main === module) {
  main();
}
