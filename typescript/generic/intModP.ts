import { IField } from './iField';
import { IMath } from './iMath';
import { IOrdered } from './iOrdered';

export class IntModP implements IField<IntModP>, IMath<IntModP>, IOrdered<IntModP> {
  value: number;
  modulus: number;
  constructor(value: number, modulus: number) {
    this.value = ((value % modulus) + modulus) % modulus;
    this.modulus = modulus;
  }
  a(o: IntModP): IntModP { return new IntModP((this.value + o.value) % this.modulus, this.modulus); }
  ae(o: IntModP): void { this.value = (this.value + o.value) % this.modulus; }
  s(o: IntModP): IntModP { return new IntModP((this.value - o.value + this.modulus) % this.modulus, this.modulus); }
  se(o: IntModP): void { this.value = (this.value - o.value + this.modulus) % this.modulus; }
  m(o: IntModP): IntModP { return new IntModP((this.value * o.value) % this.modulus, this.modulus); }
  me(o: IntModP): void { this.value = (this.value * o.value) % this.modulus; }
  d(o: IntModP): IntModP { return new IntModP((this.value * IntModP.modInv(o.value, this.modulus)) % this.modulus, this.modulus); }
  de(o: IntModP): void { this.value = (this.value * IntModP.modInv(o.value, this.modulus)) % this.modulus; }
  coerce(o: number): IntModP { return new IntModP(o, this.modulus); }
  is_zero(): boolean { return this.value === 0; }
  is_one(): boolean { return this.value === 1; }
  zero(): IntModP { return new IntModP(0, this.modulus); }
  one(): IntModP { return new IntModP(1, this.modulus); }
  abs(): number { return Math.abs(this.value); }
  sqrt(): IntModP { throw new Error('sqrt not implemented for IntModP'); }
  lt(o: IntModP): boolean { return this.value < o.value; }
  le(o: IntModP): boolean { return this.value <= o.value; }
  gt(o: IntModP): boolean { return this.value > o.value; }
  ge(o: IntModP): boolean { return this.value >= o.value; }
  e(o: IntModP): boolean { return this.value === o.value; }
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
