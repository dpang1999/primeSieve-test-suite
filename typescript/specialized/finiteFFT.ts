
import { find_prime_congruent_one_mod_n } from "../helpers/find_prime";
import { LCG } from "../helpers/lcg";
class FFT {
  static mod = 7;

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

  static transform_internal(data: number[], direction: number): void {
    if (data.length === 0) return;
    const n = data.length;
    if (n === 1) return;
    const logn = FFT.log2(n);
    FFT.bitreverse(data);
    const roots = FFT.precomputeRootsOfUnity(n, direction);
 
    let dual = 1;
    for (let bit = 0; bit < logn; bit++) {
      for (let a = 0; a < dual; a++) {
        const w = roots[a * (n / (2 * dual))];
        for (let b = 0; b < n; b += 2 * dual) {
          const i = b + a;
          const j = b + a + dual;
          const wd = w*data[j] % FFT.mod;
          const u = data[i];
          data[j] = (u + FFT.mod - wd) % FFT.mod;
          data[i] = (u + wd) % FFT.mod;
        }
      }
      dual *= 2;
    }
    
  }

  transform(data: number[]): void {
    FFT.transform_internal(data, -1);
  }

  inverse(data: number[]): void {
    FFT.transform_internal(data, 1);
    const n = data.length;
    const norm = FFT.modInv(n, FFT.mod);
    for (let d = 0; d < data.length; d++) {
      data[d] = (data[d] * norm) % FFT.mod;
    }
  }

  test(data: number[]): number {
    const nd = data.length;
    const copy = data.slice();
    this.transform(data);
    //console.log('After transform:', data);
    this.inverse(data);
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
    const nd = n;
    const data: number[] = [];
    for (let i = 0; i < nd; i++) {
      data.push(rand.nextInt() % FFT.mod);
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
  static primitive_root(): number {
    let factors = this.factorize(FFT.mod - 1);
    for (let g = 2; g < FFT.mod; g++) {
      let isPrimitiveRoot = true;
      for (let factor of factors) {
        if (this.mod_pow(g, (FFT.mod - 1) / factor, FFT.mod) === 1) {
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
  static precomputeRootsOfUnity(n: number, direction: number): number[] {
    if (FFT.mod === 0 || n <= 0) {
      throw new Error('Invalid input for precomputeRootsOfUnity');
    }
    if ((FFT.mod - 1) % n !== 0) {
      throw new Error(`Modulus ${FFT.mod} minus one must be divisible by n ${n} for precomputeRootsOfUnity`);
    }
    let root = this.primitive_root();
    let omega = this.mod_pow(root, (FFT.mod - 1) / n, FFT.mod);
    let roots: number[] = [];
    for (let k = 0; k < n; k++) {
      let exponent = (k * direction + (FFT.mod-1)) % (FFT.mod -1)
      if (exponent < 0) exponent += (FFT.mod - 1);
      roots.push(this.mod_pow(omega, exponent, FFT.mod));
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
    const test = 0; // Set to 1 to run manual test, 0 to run performance test
    if (!test) {
        const n = parseInt(process.argv[2] ?? "16", 10);
        //const prime = parseInt(process.argv[3] ?? "40961", 10);
        let prime: number;
        switch (n) {
          case 1048576:
            prime = 7340033;
            break;
          case 16777216:
            prime = 167772161;
            break;
          case 67108864:
            prime = 469762049;
            break;
          default:
            prime = find_prime_congruent_one_mod_n(n);
            break;
        }
        FFT.mod = prime;
        const fft = new FFT();
        const data = fft.make_random(n);
        console.log("Typescript Specialized finite field FFT, n=" + n)
        for (let i = 0; i < 10; i++) {
            FFT.transform_internal(data, -1);
            FFT.transform_internal(data, 1);
            console.log(`Loop ${i} done`);
        }
    }
    else {
        const in1 = [38, 0, 44, 87, 6, 45, 22, 93, 0, 0, 0, 0, 0, 0, 0, 0];
        const in2 = [80, 18, 62, 90, 17, 96, 27, 97, 0, 0, 0, 0, 0, 0, 0, 0];
        //let out = [3040, 684, 5876, 11172, 5420, 16710, 12546, 20555, 16730, 15704, 21665, 5490, 13887, 4645, 9021, 0];
        let prime = 40961;
        FFT.mod = prime;
        const data1 = [];
        const data2 = [];
        for (let i = 0; i < in1.length; i++) {
          data1.push(in1[i]);
          data2.push(in2[i]);
        }
        const fft = new FFT();
        fft.transform(data1);
        fft.transform(data2);
        
        for (let i = 0; i < data1.length; i++) {
          data1[i] = (data1[i] * data2[i]) % FFT.mod;
        }

        fft.inverse(data1);
        console.log("Result of convolution:", data1);
    }
  // print array
  //console.log(data);
  //console.log(`n=${n} => RMS Error=${rms}`);
}

if (require.main === module) {
  main();
}
