import { LCG } from "../helpers/lcg";
export type Term = {
  coefficient: number;
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

class Polynomial {
  terms: Term[];
  constructor(terms: Term[], order: TermOrder) {
    // Remove near-zero coefficients, round, and sort
    this.terms = terms
      .filter(t => Math.abs(t.coefficient) > 1e-2)
      .map(t => ({
        coefficient: Math.round(t.coefficient * 1e5) / 1e5,
        exponents: t.exponents.slice(),
      }));
    this.terms.sort((a, b) => -compareExponents(a.exponents, b.exponents, order));
  }

  add(other: Polynomial, order: TermOrder): Polynomial {
    const result = this.terms.map(t => ({ ...t }));
    for (const t of other.terms) {
      let found = false;
      for (const rt of result) {
        if (rt.exponents.length === t.exponents.length && rt.exponents.every((v, i) => v === t.exponents[i])) {
          rt.coefficient += t.coefficient;
          found = true;
          break;
        }
      }
      if (!found) result.push({ ...t });
    }
    return new Polynomial(result, order);
  }

  subtract(other: Polynomial, order: TermOrder): Polynomial {
    const result = this.terms.map(t => ({ ...t }));
    for (const t of other.terms) {
      let found = false;
      for (const rt of result) {
        if (rt.exponents.length === t.exponents.length && rt.exponents.every((v, i) => v === t.exponents[i])) {
          rt.coefficient -= t.coefficient;
          found = true;
          break;
        }
      }
      if (!found) result.push({ coefficient: -t.coefficient, exponents: t.exponents.slice() });
    }
    return new Polynomial(result, order);
  }

  multiplyByTerm(term: Term, order: TermOrder): Polynomial {
    const terms = this.terms.map(t => ({
      coefficient: t.coefficient * term.coefficient,
      exponents: t.exponents.map((v, i) => v + term.exponents[i]),
    }));
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
        if (lead.exponents.length !== divLead.exponents.length) continue;
        if (lead.exponents.every((v, i) => v >= divLead.exponents[i])) {
          const coeff = lead.coefficient / divLead.coefficient;
          const exps = lead.exponents.map((v, i) => v - divLead.exponents[i]);
          const reductionTerm: Term = { coefficient: coeff, exponents: exps };
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
    const lcmExps = lead1.exponents.map((v, i) => Math.max(v, lead2.exponents[i]));
    const scale1 = lcmExps.map((v, i) => v - lead1.exponents[i]);
    const scale2 = lcmExps.map((v, i) => v - lead2.exponents[i]);
    const scaled1 = p1.multiplyByTerm({ coefficient: 1, exponents: scale1 }, order);
    const scaled2 = p2.multiplyByTerm({ coefficient: 1, exponents: scale2 }, order);
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
    const reduced = poly.reduce(basisExcl, order);
    const key = JSON.stringify(reduced.terms);
    if (reduced.terms.length > 0 && !reducedBasis.some(rb => JSON.stringify(rb.terms) === key)) {
      reducedBasis.push(reduced);
    }
  }
  return reducedBasis;
}


function main() {
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
        const coeff = rand.nextDouble() * 2 - 1;
        const exps = [rand.nextInt() % 4, rand.nextInt() % 4, rand.nextInt() % 4];
        terms.push({ coefficient: coeff, exponents: exps });
      }
      polys.push(new Polynomial(terms, order));
    }
    const basis = naiveGrobnerBasis(polys, order);
    console.log(basis.length);
    //console.log('Computed Grobner Basis:', basis);
    return;
  } else {
    const order = TermOrder.Lex;
    // Example: x^3 + y^3 + z^3
    const p1 = new Polynomial([
      { coefficient: 1, exponents: [3, 0, 0] },
      { coefficient: 1, exponents: [0, 3, 0] },
      { coefficient: 1, exponents: [0, 0, 3] },
    ], order);
    // xy + yz + xz
    const p2 = new Polynomial([
      { coefficient: 1, exponents: [1, 1, 0] },
      { coefficient: 1, exponents: [0, 1, 1] },
      { coefficient: 1, exponents: [1, 0, 1] },
    ], order);
    // x+y+z
    const p3 = new Polynomial([
      { coefficient: 1, exponents: [1, 0, 0] },
      { coefficient: 1, exponents: [0, 1, 0] },
      { coefficient: 1, exponents: [0, 0, 1] },
    ], order);
    const basis = naiveGrobnerBasis([p1, p2, p3], order);
    console.log('Final Grobner Basis:');
    for (const poly of basis) {
      for (const term of poly.terms) {
        process.stdout.write(`${term.coefficient}*[${term.exponents[0]}, ${term.exponents[1]}, ${term.exponents[2]}] + `);
      }
      console.log();
    }
  }
}

if (require.main === module) {
  main();
}