// VecExponents: vector-based exponents for generic polynomials
export class VecExponents {
  exponents: number[];
  constructor(exponents: number[]) {
    this.exponents = [...exponents];
  }
  add(o: VecExponents): VecExponents {
    return new VecExponents(this.exponents.map((v, i) => v + o.exponents[i]));
  }
  sub(o: VecExponents): VecExponents {
    return new VecExponents(this.exponents.map((v, i) => v - o.exponents[i]));
  }
  lcm(o: VecExponents): VecExponents {
    return new VecExponents(this.exponents.map((v, i) => Math.max(v, o.exponents[i])));
  }
  deg(): number {
    return this.exponents.reduce((a, b) => a + b, 0);
  }
  lexCompare(o: VecExponents): number {
    for (let i = 0; i < this.exponents.length; i++) {
      if (this.exponents[i] < o.exponents[i]) return -1;
      if (this.exponents[i] > o.exponents[i]) return 1;
    }
    return 0;
  }
  canReduce(o: VecExponents): boolean {
    return this.exponents.every((v, i) => v >= o.exponents[i]);
  }
  equals(o: VecExponents): boolean {
    return this.exponents.length === o.exponents.length && this.exponents.every((v, i) => v === o.exponents[i]);
  }
  toString(): string {
    return `[${this.exponents.join(", ")}]`;
  }
}
