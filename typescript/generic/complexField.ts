import { IField } from './iField';
import { IMath } from './iMath';
import { IOrdered } from './iOrdered';
import { IPrimitiveRoots } from './iPrimitiveRoots';

// Generic ComplexField<T> where T is a field (e.g., DoubleField, IntModP, etc.)
export class ComplexField<T extends IField<T> & IMath<T> & IOrdered<T>> implements IField<ComplexField<T>>, IMath<ComplexField<T>>, IOrdered<ComplexField<T>>, IPrimitiveRoots<ComplexField<T>>  {
  re: T;
  im: T;

  constructor(re: T, im: T) {
    this.re = re;
    this.im = im;
  }

  a(o: ComplexField<T>): ComplexField<T> {
    return new ComplexField(this.re.a(o.re), this.im.a(o.im));
  }
  ae(o: ComplexField<T>): void {
    this.re.ae(o.re);
    this.im.ae(o.im);
  }
  s(o: ComplexField<T>): ComplexField<T> {
    return new ComplexField(this.re.s(o.re), this.im.s(o.im));
  }
  se(o: ComplexField<T>): void {
    this.re.se(o.re);
    this.im.se(o.im);
  }
  m(o: ComplexField<T>): ComplexField<T> {
    // (a + bi) * (c + di) = (ac - bd) + (ad + bc)i
    const real = this.re.m(o.re).s(this.im.m(o.im));
    const imag = this.re.m(o.im).a(this.im.m(o.re));
    return new ComplexField(real, imag);
  }
  me(o: ComplexField<T>): void {
    const real = this.re.m(o.re).s(this.im.m(o.im));
    const imag = this.re.m(o.im).a(this.im.m(o.re));
    this.re = real;
    this.im = imag;
  }
  d(o: ComplexField<T>): ComplexField<T> {
    // (a + bi) / (c + di) = [(ac + bd) / (c^2 + d^2)] + [(bc - ad) / (c^2 + d^2)]i
    const denom = o.re.m(o.re).a(o.im.m(o.im));
    const real = this.re.m(o.re).a(this.im.m(o.im)).d(denom);
    const imag = this.im.m(o.re).s(this.re.m(o.im)).d(denom);
    return new ComplexField(real, imag);
  }
  de(o: ComplexField<T>): void {
    const denom = o.re.m(o.re).a(o.im.m(o.im));
    const real = this.re.m(o.re).a(this.im.m(o.im)).d(denom);
    const imag = this.im.m(o.re).s(this.re.m(o.im)).d(denom);
    this.re = real;
    this.im = imag;
  }
  coerce(o: number): ComplexField<T> {
    return new ComplexField(this.re.coerce(o), this.im.zero());
  }
  coerce_to_number(): number {
    return Math.pow(this.re.coerce_to_number(), 2) + Math.pow(this.im.coerce_to_number(), 2);
  }
  is_zero(): boolean {
    return this.re.is_zero() && this.im.is_zero();
  }
  is_one(): boolean {
    return this.re.is_one() && this.im.is_zero();
  }
  zero(): ComplexField<T> {
    return new ComplexField(this.re.zero(), this.im.zero());
  }
  one(): ComplexField<T> {
    return new ComplexField(this.re.one(), this.im.zero());
  }
  abs(): number{
    // sqrt(re^2 + im^2)
    return Math.sqrt(this.re.m(this.re).a(this.im.m(this.im)).abs());
  }
  sqrt(): void {
    // sqrt(a + bi) = sqrt((|z| + a)/2) + sign(b) * sqrt((|z| - a)/2) i
    const modulus = this.re.m(this.re).a(this.im.m(this.im));
    const two = this.re.coerce(2);
    const real = modulus.a(this.re).d(two);
    real.sqrt();
    let imag = modulus.s(this.re).d(two);
    imag.sqrt();
    if (this.im.lt(this.im.zero())) {
      imag = imag.m(this.re.coerce(-1));
    }
    this.re = real;
    this.im = imag;
  }
  lt(o: ComplexField<T>): boolean {
    // Compare real, then imaginary
    if (this.re.lt(o.re)) return true;
    if (this.re.eq(o.re)) return this.im.lt(o.im);
    return false;
  }
  le(o: ComplexField<T>): boolean {
    if (this.re.lt(o.re)) return true;
    if (this.re.eq(o.re)) return this.im.le(o.im);
    return false;
  }
  gt(o: ComplexField<T>): boolean {
    if (this.re.gt(o.re)) return true;
    if (this.re.eq(o.re)) return this.im.gt(o.im);
    return false;
  }
  ge(o: ComplexField<T>): boolean {
    if (this.re.gt(o.re)) return true;
    if (this.re.eq(o.re)) return this.im.ge(o.im);
    return false;
  }
  eq(o: ComplexField<T>): boolean {
    return this.re.eq(o.re) && this.im.eq(o.im);
  }
  copy(): ComplexField<T> {
    return new ComplexField(this.re.copy(), this.im.copy());
  }
  toString(): string {
    const reStr = this.re.toString();
    const imStr = this.im.toString();
    if (this.im.is_zero()) return `(${reStr})`;
    if (this.im.lt(this.im.zero())) return `(${reStr}${imStr}i)`;
    return `(${reStr}+${imStr}i)`;
  }


  primitive_root(n: number): ComplexField<T> {
    if (n <= 0) {
      throw new Error("n must be positive");
    }

    // If T supports trigonometric functions (duck-typing)
    const trigRe = this.re as any;
    if (typeof trigRe.coerce === 'function' && typeof trigRe.cos === 'function' && typeof trigRe.sin === 'function') {
    // Compute the angle for the primitive root of unity
    const angle = trigRe.coerce(2 * Math.PI / n);
    const realPart = angle.cos();
    const imagPart = angle.sin();
    return new ComplexField(realPart, imagPart);
    } else if (typeof trigRe.primitive_root === 'function') {
    // For finite fields, compute the primitive root algebraically
    const realPart = trigRe.primitive_root(n);
    const imagPart = this.im.zero();
    return new ComplexField(realPart, imagPart);
    }
    throw new Error("Unsupported field type for primitive root calculation");
}


  pow(exp: number): ComplexField<T> {
    if (exp === 0) return this.one();
    if (exp < 0) return this.inverse().pow(-exp);

    // Convert to polar form
    const r = this.coerce_to_number();

    // Actually, we want the numeric value of re and im for atan2:
    const theta2 = Math.atan2(this.im.coerce_to_number(), this.re.coerce_to_number());
    const newR = Math.pow(r, exp);
    const newTheta = theta2 * exp;

    const real = this.re.coerce(newR * Math.cos(newTheta));
    const imag = this.im.coerce(newR * Math.sin(newTheta));
    return new ComplexField(real, imag);
  }

  precomputeRootsOfUnity(n: number, direction: number): ComplexField<T>[] {
    const roots: ComplexField<T>[] = [];
    for (let k = 0; k < n; k++) {
      const angle = 2.0 * Math.PI * k * direction / n;
      const realPart = this.re.coerce(Math.cos(angle));
      const imagPart = this.im.coerce(Math.sin(angle));
      roots.push(new ComplexField(realPart, imagPart));
    }
    return roots;
  }

  inverse(): ComplexField<T> {
    const denom = this.re.m(this.re).a(this.im.m(this.im));
    const real = this.re.d(denom);
    const imag = this.im.m(this.re.coerce(-1)).d(denom);
    return new ComplexField(real, imag);
  }
}