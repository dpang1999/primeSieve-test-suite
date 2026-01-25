
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
    const n = data.length / 2;
    const nm1 = n - 1;
    let i = 0, j = 0;
    while (i < nm1) {
      if (i < j) {
        const ii = i * 2;
        const jj = j * 2;
        [data[ii], data[jj]] = [data[jj], data[ii]];
        [data[ii + 1], data[jj + 1]] = [data[jj + 1], data[ii + 1]];
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
    const n = data.length / 2;
    if (n === 1) return;
    const logn = FFT.log2(n);
    FFT.bitreverse(data);
    for (let bit = 0; bit < logn; bit++) {
      const dual = 1 << bit;
      let w_real = 1.0;
      let w_imag = 0.0;
      const theta = 2.0 * direction * Math.PI / (2.0 * dual);
      const s = Math.sin(theta);
      const t = Math.sin(theta / 2.0);
      const s2 = 2.0 * t * t;
      // a = 0
      for (let b = 0; b < n; b += 2 * dual) {
        const i = 2 * b;
        const j = 2 * (b + dual);
        const wd_real = data[j];
        const wd_imag = data[j + 1];
        data[j] = data[i] - wd_real;
        data[j + 1] = data[i + 1] - wd_imag;
        data[i] += wd_real;
        data[i + 1] += wd_imag;
      }
      // a = 1 .. (dual-1)
      for (let a = 1; a < dual; a++) {
        const tmp_real = w_real - s * w_imag - s2 * w_real;
        const tmp_imag = w_imag + s * w_real - s2 * w_imag;
        w_real = tmp_real;
        w_imag = tmp_imag;
        for (let b = 0; b < n; b += 2 * dual) {
          const i = 2 * (b + a);
          const j = 2 * (b + a + dual);
          const z1_real = data[j];
          const z1_imag = data[j + 1];
          const wd_real = w_real * z1_real - w_imag * z1_imag;
          const wd_imag = w_real * z1_imag + w_imag * z1_real;
          data[j] = data[i] - wd_real;
          data[j + 1] = data[i + 1] - wd_imag;
          data[i] += wd_real;
          data[i + 1] += wd_imag;
        }
      }
    }
  }

  transform(data: number[]): void {
    FFT.transform_internal(data, -1);
  }

  inverse(data: number[]): void {
    FFT.transform_internal(data, 1);
    const n = data.length / 2;
    const norm = 1.0 / n;
    for (let d = 0; d < data.length; d++) {
      data[d] *= norm;
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
    const nd = 2 * n;
    const data: number[] = [];
    for (let i = 0; i < nd; i++) {
      data.push(rand.nextDouble());
    }
    return data;
  }
}


function main() {
  const n = parseInt(process.argv[2] ?? "16", 10);
  const fft = new FFT();
  const data = fft.make_random(n);
  const rms = fft.test(data);
  // print array
  //console.log(data);
  console.log(`n=${n} => RMS Error=${rms}`);
}

if (require.main === module) {
  main();
}
