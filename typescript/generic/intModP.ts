import { IField } from './iField';
import { IMath } from './iMath';
import { IOrdered } from './iOrdered';
import { IPrimitiveRoots } from './iPrimitiveRoots';

function modPow(base: number, exp: number, modulus: number): number {
  let result = 1;
  base = base % modulus;
  while (exp > 0) {
    if (exp % 2 === 1) result = (result * base) % modulus;
    base = (base * base) % modulus;
    exp = Math.floor(exp / 2);
  }
  return result;
}

function factorize(n: number): number[] {
  const factors: number[] = [];
  let i = 2;
  while (i * i <= n) {
    if (n % i === 0) {
      factors.push(i);
      while (n % i === 0) n = Math.floor(n / i);
    }
    i++;
  }
  if (n > 1) factors.push(n);
  return factors;
}

export class IntModP implements IField<IntModP>, IMath<IntModP>, IOrdered<IntModP>, IPrimitiveRoots<IntModP> {
    // Returns a primitive root modulo p (for p prime)
    primitive_root(n: number): IntModP {
      if (n === 0 || n >= this.modulus) {
        throw new Error('n must be in range [1, p-1]');
      }
      const p = this.modulus;
      const factors = factorize(p - 1);
      for (let g = 2; g < p; g++) {
        let isRoot = true;
        for (const factor of factors) {
          if (modPow(g, (p - 1) / factor, p) === 1) {
            isRoot = false;
            break;
          }
        }
        if (isRoot) {
          return new IntModP(g, p);
        }
      }
      return new IntModP(0, p); // No primitive root found
    }

    pow(exp: number): IntModP {
      return new IntModP(modPow(this.value, exp, this.modulus), this.modulus);
    }

    precomputeRootsOfUnity(n: number, direction: number): IntModP[] {
      if ((this.modulus - 1) % n !== 0) {
        throw new Error('n must divide p-1 for roots of unity to exist in IntModP');
      }
      const p = this.modulus;
      const g = this.primitive_root(p - 1);
      const omega = g.pow((p - 1) / n);
      const roots: IntModP[] = [];
      for (let k = 0; k < n; k++) {
        let exponent = (k * direction) % (p - 1);
        if (exponent < 0) exponent += (p - 1);
        roots.push(omega.pow(exponent));
      }
      return roots;
    }
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
  coerce_to_number(): number { return this.value; }
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
