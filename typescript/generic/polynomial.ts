import { IField } from './iField';
import { IMath } from './iMath';
import { IOrdered } from './iOrdered';
import { ICopiable } from './iCopiable';

// Generic Term and Polynomial for any field and exponent type
export class Term<C, E> {
  coefficient: C;
  exponents: E;
  constructor(coefficient: C, exponents: E) {
    this.coefficient = coefficient;
    this.exponents = exponents;
  }
  toString(): string {
    return `${this.coefficient}*${this.exponents}`;
  }
}

export class Polynomial<C, E> {
  terms: Term<C, E>[];
  constructor(terms: Term<C, E>[]) {
    this.terms = terms;
  }
  toString(): string {
    if (this.terms.length === 0) return '0';
    return this.terms.map(t => t.toString()).join(' + ');
  }
}
