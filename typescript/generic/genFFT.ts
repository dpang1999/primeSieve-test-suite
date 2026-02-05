import { find_prime_congruent_one_mod_n } from '../helpers/find_prime';
import { LCG } from '../helpers/lcg';
import { IField } from './iField';
import { IMath } from './iMath';
import { IOrdered } from './iOrdered';
import { IPrimitiveRoots } from './iPrimitiveRoots';

// Class-based generic FFT, similar to the Rust GenFFT
export class GenFFT<C extends IField<C> & IOrdered<C> & IPrimitiveRoots<C> & IMath<C>> {
  c: C;

  constructor(data: C) {
    this.c = data;
  }

  transform(data: C[]): void {
    this.transformInternal(data, -1);
  }

  inverse(data: C[]): void {
    this.transformInternal(data, 1);
    const n = data.length;
    const norm = this.c.coerce(n);
    for (let i = 0; i < n; i++) {
      data[i].de(norm);
    }
  }

  test(data: C[]): number {
    const n = data.length;
    const copy = data.map(x => x.copy());

    //print before transform
    console.log("Before transform:", data.map(d => d.toString()).join(", "));

    this.transform(data);
    //print after transform
    console.log("After transform:", data.map(d => d.toString()).join(", "));

    this.inverse(data);
    //print after inverse
    console.log("After inverse:", data.map(d => d.toString()).join(", "));

    let diff = 0.0;
    for (let i = 0; i < n; i++) {
      const d = data[i];
      const orig = copy[i];
      const realDiff = d.coerce_to_number() - orig.coerce_to_number();
      diff += realDiff * realDiff;
    }
    return Math.sqrt(diff / n);
  }

  private transformInternal(data: C[], direction: number): void {
    const n = data.length;
    if (n === 0 || n === 1) return;

    this.bitreverse(data);

    const roots = this.c.precomputeRootsOfUnity(n, direction);

    let dual = 1;
    const logn = Math.log2(n);
    for (let bit = 0; bit < logn; bit++) {
      for (let a = 0; a < dual; a++) {
        const w = roots[a * (n / (2 * dual))];
        for (let b = 0; b < n; b += 2 * dual) {
          const i = b + a;
          const j = b + a + dual;
          const wd = w.m(data[j]);
          const u = data[i];
          data[j] = u.s(wd);
          data[i] = u.a(wd);
        }
      }
      dual *= 2;
    }
  }

  private bitreverse(data: C[]): void {
    const n = data.length;
    let j = 0;
    for (let i = 0; i < n; i++) {
      if (i < j) [data[i], data[j]] = [data[j], data[i]];
      let k = n >> 1;
      while (k && j >= k) {
        j -= k;
        k >>= 1;
      }
      j += k;
    }
  }
}

if(require.main == module) {
  // Usage: tsx genFFT.ts [n] [mode]
  // arg1 = size (N = power of 2)
  // arg2 = field type (0 = finite field, 1 = complex field)
  const args = process.argv.slice(2);
  const n = parseInt(args[0] || '16', 10);
  const mode = parseInt(args[1] || '0', 10);
  const rand = new LCG(12345, 1345, 65, 17);

  if (mode === 1) {
    const { ComplexField } = require('./complexField');
    const { DoubleField } = require('./doubleField');
    const data = [];
    for (let i = 0; i < n; i++) {
      data.push(new ComplexField(new DoubleField(rand.nextDouble()), new DoubleField(rand.nextDouble())));
    }
    const fft = new GenFFT(data[0]);
    const error = fft.test(data);

    console.log(`FFT test completed with RMS error: ${error}`);
  } else {
    const { IntModP } = require('./intModP');
    const data = [];
    const prime = find_prime_congruent_one_mod_n(n);
    IntModP.setModulus(prime);
    for (let i = 0; i < n; i++) {
      data.push(new IntModP(rand.nextInt() % prime));
    }
    const fft = new GenFFT(data[0]);
    const error = fft.test(data);
    console.log(`FFT test completed with RMS error: ${error}`);
  } 
}
