// BitPackedExponents: bit-packed exponents for up to 6 variables
export class BitPackedExponents {
  packed: bigint;
  constructor(packed: bigint) {
    this.packed = packed;
  }
  static fromArray(exps: number[]): BitPackedExponents {
    let packed = 0n;
    for (let i = 0; i < 6; i++) {
      const shift = BigInt(40 - 8 * i);
      packed |= BigInt(exps[i] || 0) << shift;
    }
    const degree = exps.reduce((a, b) => a + (b || 0), 0);
    packed |= BigInt(degree) << 48n;
    return new BitPackedExponents(packed);
  }
  toArray(): number[] {
    const arr = [];
    for (let i = 0; i < 6; i++) {
      const shift = BigInt(40 - 8 * i);
      arr.push(Number((this.packed >> shift) & 0xFFn));
    }
    return arr;
  }
  deg(): number {
    return Number((this.packed >> 48n) & 0xFFFFn);
  }
  add(o: BitPackedExponents): BitPackedExponents {
    return new BitPackedExponents(this.packed + o.packed);
  }
  sub(o: BitPackedExponents): BitPackedExponents {
    return new BitPackedExponents(this.packed - o.packed);
  }
  lcm(o: BitPackedExponents): BitPackedExponents {
    let lcm = 0n;
    let degree = 0;
    for (let i = 0; i < 6; i++) {
      const shift = BigInt(40 - 8 * i);
      const a = (this.packed >> shift) & 0xFFn;
      const b = (o.packed >> shift) & 0xFFn;
      const l = a > b ? a : b;
      lcm |= l << shift;
      degree += Number(l);
    }
    lcm |= BigInt(degree) << 48n;
    return new BitPackedExponents(lcm);
  }
  lexCompare(o: BitPackedExponents): number {
    for (let i = 0; i < 6; i++) {
      const shift = BigInt(40 - 8 * i);
      const a = (this.packed >> shift) & 0xFFn;
      const b = (o.packed >> shift) & 0xFFn;
      if (a < b) return -1;
      if (a > b) return 1;
    }
    return 0;
  }
  canReduce(o: BitPackedExponents): boolean {
    for (let i = 0; i < 6; i++) {
      const shift = BigInt(40 - 8 * i);
      const a = (this.packed >> shift) & 0xFFn;
      const b = (o.packed >> shift) & 0xFFn;
      if (a < b) return false;
    }
    return true;
  }
  equals(o: BitPackedExponents): boolean {
    return this.packed === o.packed;
  }
  toString(): string {
    return `0x${this.packed.toString(16)}`;
  }
}
