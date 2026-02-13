import { IField } from './iField';
import { IMath } from './iMath';
import { IOrdered } from './iOrdered';

export class IntModP implements IField<IntModP>, IMath<IntModP>, IOrdered<IntModP> {
  value: number;
  static modulus: number;
  static setModulus(modulus: number) {
    if (modulus <= 1) throw new Error('Modulus must be greater than 1');
    IntModP.modulus = modulus;
  }
  static getModulus(): number {
    if (IntModP.modulus === undefined) throw new Error('Modulus not set');
    return IntModP.modulus;
  }
  constructor(value: number) {
    this.value = ((value % IntModP.modulus) + IntModP.modulus) % IntModP.modulus;
  }
  a(o: IntModP): IntModP { return new IntModP((this.value + o.value) % IntModP.modulus); }
  ae(o: IntModP): void { this.value = (this.value + o.value) % IntModP.modulus; }
  s(o: IntModP): IntModP { return new IntModP((this.value - o.value + IntModP.modulus) % IntModP.modulus); }
  se(o: IntModP): void { this.value = (this.value - o.value + IntModP.modulus) % IntModP.modulus; }
  m(o: IntModP): IntModP { return new IntModP((this.value * o.value) % IntModP.modulus); }
  me(o: IntModP): void { this.value = (this.value * o.value) % IntModP.modulus; }
  d(o: IntModP): IntModP { return new IntModP((this.value * IntModP.modInv(o.value, IntModP.modulus)) % IntModP.modulus); }
  de(o: IntModP): void { this.value = (this.value * IntModP.modInv(o.value, IntModP.modulus)) % IntModP.modulus; }
  coerce(o: number): IntModP { return new IntModP(o); }
  coerce_to_number(): number { return this.value; }
  is_zero(): boolean { return this.value === 0; }
  is_one(): boolean { return this.value === 1; }
  zero(): IntModP { return new IntModP(0); }
  one(): IntModP { return new IntModP(1); }
  abs(): number { return Math.abs(this.value); }
  sqrt(): IntModP { throw new Error('sqrt not implemented for IntModP'); }
  lt(o: IntModP): boolean { return this.value < o.value; }
  le(o: IntModP): boolean { return this.value <= o.value; }
  gt(o: IntModP): boolean { return this.value > o.value; }
  ge(o: IntModP): boolean { return this.value >= o.value; }
  eq(o: IntModP): boolean { return this.value === o.value; }
  copy(): IntModP { return new IntModP(this.value); }
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
  primitive_root(p: number): IntModP {
    let mod = IntModP.getModulus();
    let factors = this.factorize(p-1);
    for (let g = 2; g < mod; g++) {
      let isPrimitiveRoot = true;
      for (let factor of factors) {
        if (this.mod_pow(g, (mod - 1) / factor) === 1) {
          isPrimitiveRoot = false;
          break;
        }
      }
      if (isPrimitiveRoot) {
        return new IntModP(g);
      }
    }
    throw new Error('No primitive root found');
  }
  pow(exp: number): IntModP {
    return new IntModP(this.mod_pow(this.value, exp));
  }
  precomputeRootsOfUnity(n: number, direction: number): IntModP[] {
    let mod = IntModP.getModulus();
    if ((mod - 1) % n !== 0) {
      throw new Error('Modulus minus one must be divisible by n for precomputeRootsOfUnity');
    }
    let root = this.primitive_root(mod);
    let omega = root.pow((mod - 1) / n);
    let roots: IntModP[] = [];
    for (let k = 0; k < n; k++) {
      let exponent = (k * direction + (mod-1)) % (mod -1)
      if (exponent < 0) exponent += (mod - 1);
      roots.push(omega.pow(exponent));
    }
    return roots;
  }
  factorize(n: number): number[] {
    let factors: number[] = [];
    for (let i = 2; i * i <= n; i++) {
      if (n % i === 0) {
        factors.push(i);
        while (n % i === 0) {
          n = Math.floor(n / i);
        }
      }
    }
    if (n > 1) {
      factors.push(n);
    }
    return factors;
  }
  mod_pow(base: number, exp: number): number {
    let result = 1;
    let modulus = IntModP.getModulus();
    while (exp > 0) {
      if (exp % 2 === 1) {
        result = (result * base) % modulus;
      }
      exp = Math.floor(exp / 2);
      base = (base * base) % modulus;
    }
    return result;
  }
  
}
