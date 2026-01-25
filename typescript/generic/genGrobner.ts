import { IField } from './iField';
import { IMath } from './iMath';
import { IOrdered } from './iOrdered';
import { ICopiable } from './iCopiable';
import { Term, Polynomial } from './polynomial';
import { TermOrder, naiveGrobnerBasis } from './grobner';
import { DoubleField } from './doubleField';
import { SingleField } from './singleField';
import { IntModP } from './intModP';
import { VecExponents } from './vecExponents';
import { BitPackedExponents } from './bitPackedExponents';

// Utility to generate a random polynomial basis
type FieldType = 'double' | 'single' | 'intmodp';
type ExponentType = 'vec' | 'bitpacked';

export function randomBasis(
  nPolys: number,
  nVars: number,
  deg: number,
  fieldType: FieldType,
  exponentType: ExponentType,
  p?: number // for IntModP
): Polynomial<any, any>[] {
  let field: IField<any>;
  switch (fieldType) {
    case 'double':
      field = new DoubleField();
      break;
    case 'single':
      field = new SingleField();
      break;
    case 'intmodp':
      if (typeof p !== 'number') throw new Error('p required for IntModP');
      field = new IntModP(p);
      break;
    default:
      throw new Error('Unknown field type');
  }
  const polys: Polynomial<any, any>[] = [];
  for (let i = 0; i < nPolys; i++) {
    const terms: Term<any, any>[] = [];
    const nTerms = Math.floor(Math.random() * nVars) + 1;
    for (let j = 0; j < nTerms; j++) {
      let coeff = field.random();
      let exponents;
      switch (exponentType) {
        case 'vec':
          exponents = VecExponents.random(nVars, deg);
          break;
        case 'bitpacked':
          exponents = BitPackedExponents.random(nVars, deg);
          break;
        default:
          throw new Error('Unknown exponent type');
      }
      terms.push(new Term(coeff, exponents));
    }
    polys.push(new Polynomial(terms));
  }
  return polys;
}

// Example test harness
type GrobnerTestArgs = {
  nPolys: number;
  nVars: number;
  deg: number;
  fieldType: FieldType;
  exponentType: ExponentType;
  p?: number;
  order?: TermOrder;
};

export function runGrobnerTest(args: GrobnerTestArgs) {
  const order = args.order ?? TermOrder.Lex;
  const basis = randomBasis(
    args.nPolys,
    args.nVars,
    args.deg,
    args.fieldType,
    args.exponentType,
    args.p
  );
  const grobner = naiveGrobnerBasis(basis, order);
  console.log('Computed Grobner Basis Polynomials:');
  for (const poly of grobner) {
    console.log(poly.toString());
  }
}
