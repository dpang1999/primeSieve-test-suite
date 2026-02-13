
import { LCG } from "../helpers/lcg";
class FFT {
  static log2(n: number): number {
    let log = 0;
    let k = 1;
    while (k < n) {
      k *= 2;
      log += 1;
    }
    if (n !== (1 << log)) throw new Error(`FFT: Data length is not a power of 2!: ${n}`);
    return log;
  }

  static bitreverse(data: number[]): void {
    const n = data.length;
    const nm1 = n - 1;
    let i = 0, j = 0;
    while (i < nm1) {
      if (i < j) {
        [data[i], data[j]] = [data[j], data[i]];
        }
      let k = n >> 1;
      while (k <= j) {
        j -= k;
        k >>= 1;
      }
      j += k;
      i += 1;
    }
  }

  static transform_internal(data: number[], direction: number, mod: number): void {
    if (data.length === 0) return;
    const n = data.length;
    if (n === 1) return;
    const logn = FFT.log2(n);
    FFT.bitreverse(data);
    const roots = FFT.precomputeRootsOfUnity(n, direction, 40961);
 
    let dual = 1;
    for (let bit = 0; bit < logn; bit++) {
      for (let a = 0; a < dual; a++) {
        const w = roots[a * (n / (2 * dual))];
        for (let b = 0; b < n; b += 2 * dual) {
          const i = b + a;
          const j = b + a + dual;
          const wd = w*data[j] % mod;
          const u = data[i];
          data[j] = (u + mod - wd) % mod;
          data[i] = (u + wd) % mod;
        }
      }
      dual *= 2;
    }
    
  }

  transform(data: number[], mod: number): void {
    FFT.transform_internal(data, -1, mod);
  }

  inverse(data: number[], mod: number): void {
    FFT.transform_internal(data, 1, mod);
    const n = data.length;
    const norm = FFT.modInv(n, mod);
    for (let d = 0; d < data.length; d++) {
      data[d] = (data[d] * norm) % mod;
    }
  }

  test(data: number[], mod: number): number {
    const nd = data.length;
    const copy = data.slice();
    this.transform(data, mod);
    //console.log('After transform:', data);
    this.inverse(data, mod);
    //console.log('After inverse:', data);
    let diff = 0.0;
    for (let i = 0; i < nd; i++) {
      const d = data[i] - copy[i];
      diff += d * d;
    }
    return Math.sqrt(diff / nd);
  }

  make_random(n: number): number[] {
    // Interleaved real/imag, like Rust
    const rand = new LCG(12345, 1345, 16645, 1013904);
    const nd = 2 * n;
    const data: number[] = [];
    for (let i = 0; i < nd; i++) {
      data.push(rand.nextDouble());
    }
    return data;
  }
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
  static primitive_root(mod: number): number {
    let factors = this.factorize(mod - 1);
    for (let g = 2; g < mod; g++) {
      let isPrimitiveRoot = true;
      for (let factor of factors) {
        if (this.mod_pow(g, (mod - 1) / factor, mod) === 1) {
          isPrimitiveRoot = false;
          break;
        }
      }
      if (isPrimitiveRoot) {
        return g
      }
    }
    throw new Error('No primitive root found');
  }
  static precomputeRootsOfUnity(n: number, direction: number, mod: number): number[] {
    if (mod === 0 || n <= 0) {
      throw new Error('Invalid input for precomputeRootsOfUnity');
    }
    if ((mod - 1) % n !== 0) {
      throw new Error('Modulus minus one must be divisible by n for precomputeRootsOfUnity');
    }
    let root = this.primitive_root(mod);
    let omega = this.mod_pow(root, (mod - 1) / n, mod);
    let roots: number[] = [];
    for (let k = 0; k < n; k++) {
      let exponent = (k * direction + (mod-1)) % (mod -1)
      if (exponent < 0) exponent += (mod - 1);
      roots.push(this.mod_pow(omega, exponent, mod));
    }
    return roots;
  }
  static factorize(n: number): number[] {
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
  static mod_pow(base: number, exp: number, mod: number): number {
    let result = 1;
    while (exp > 0) {
      if (exp % 2 === 1) {
        result = (result * base) % mod;
      }
      exp = Math.floor(exp / 2);
      base = (base * base) % mod;
    }
    return result;
  }
  
}


function main() {
    const test = 1; // Set to 1 to run the test, 0 to just run the FFT without checking error
    if (!test) {
        const n = parseInt(process.argv[2] ?? "16", 10);
        const prime = parseInt(process.argv[3] ?? "40961", 10);
        const fft = new FFT();
        const data = fft.make_random(n);
        console.log("Typescript Specialized number FFT, n=" + n)
        for (let i = 0; i < 10; i++) {
            FFT.transform_internal(data, -1, prime);
            FFT.transform_internal(data, 1, prime);
            console.log(`Loop ${i} done`);
        }
    }
    else {
        const in1 = [38, 0, 44, 87, 6, 45, 22, 93, 0, 0, 0, 0, 0, 0, 0, 0];
        const in2 = [80, 18, 62, 90, 17, 96, 27, 97, 0, 0, 0, 0, 0, 0, 0, 0];
        //let out = [3040, 684, 5876, 11172, 5420, 16710, 12546, 20555, 16730, 15704, 21665, 5490, 13887, 4645, 9021, 0];
        let prime = 40961;
        const data1 = [];
        const data2 = [];
        for (let i = 0; i < in1.length; i++) {
          data1.push(in1[i]);
          data2.push(in2[i]);
        }
        const fft = new FFT();
        fft.transform(data1, prime);
        fft.transform(data2, prime);
        
        for (let i = 0; i < data1.length; i++) {
          data1[i] = (data1[i] * data2[i]) % prime;
        }

        fft.inverse(data1, prime);
        console.log("Result of convolution:", data1);
    }
  // print array
  //console.log(data);
  //console.log(`n=${n} => RMS Error=${rms}`);
}

if (require.main === module) {
  main();
}
