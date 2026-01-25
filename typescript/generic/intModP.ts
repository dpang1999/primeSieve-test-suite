import { IField } from './iField';
import { IMath } from './iMath';
import { IOrdered } from './iOrdered';
import { ICopiable } from './iCopiable';

export class IntModP implements IField<IntModP>, IMath<IntModP>, IOrdered<IntModP>, ICopiable<IntModP> {
  value: number;
  modulus: number;
  constructor(value: number, modulus: number) {
    this.value = ((value % modulus) + modulus) % modulus;
    this.modulus = modulus;
  }
  add(o: IntModP): IntModP { return new IntModP((this.value + o.value) % this.modulus, this.modulus); }
  sub(o: IntModP): IntModP { return new IntModP((this.value - o.value + this.modulus) % this.modulus, this.modulus); }
  mul(o: IntModP): IntModP { return new IntModP((this.value * o.value) % this.modulus, this.modulus); }
  div(o: IntModP): IntModP { return new IntModP((this.value * IntModP.modInv(o.value, this.modulus)) % this.modulus, this.modulus); }
  coerceFromInt(i: number): IntModP { return new IntModP(i, this.modulus); }
  coerceFromFloat(f: number): IntModP { return new IntModP(Math.round(f), this.modulus); }
  coerceToFloat(): number { return this.value; }
  isZero(): boolean { return this.value === 0; }
  isOne(): boolean { return this.value === 1; }
  zero(): IntModP { return new IntModP(0, this.modulus); }
  one(): IntModP { return new IntModP(1, this.modulus); }
  abs(): IntModP { return new IntModP(Math.abs(this.value), this.modulus); }
  sqrt(): IntModP { throw new Error('sqrt not implemented for IntModP'); }
  lt(o: IntModP): boolean { return this.value < o.value; }
  le(o: IntModP): boolean { return this.value <= o.value; }
  gt(o: IntModP): boolean { return this.value > o.value; }
  ge(o: IntModP): boolean { return this.value >= o.value; }
  eq(o: IntModP): boolean { return this.value === o.value; }
  copy(): IntModP { return new IntModP(this.value, this.modulus); }
  toString(): string { return this.value.toString(); }
  static modInv(a: number, m: number): number {
    let m0 = m, t, q;
    let x0 = 0, x1 = 1;
    if (m === 1) return 0;
    while (a > 1) {
      q = Math.floor(a / m);
      t = m;
      m = a % m; a = t;
      t = x0;
      x0 = x1 - q * x0;
      x1 = t;
    }
    if (x1 < 0) x1 += m0;
    return x1;
  }
}
