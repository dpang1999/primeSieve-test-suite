import { IField } from './iField';
import { IMath } from './iMath';
import { IPrimitiveRoots } from './iPrimitiveRoots';

// Generic Cooley-Tukey FFT for any field with primitive roots of unity
export function genericFFT<C>(
  input: C[],
  field: IField<C> & IMath<C> & IPrimitiveRoots<C>,
  inverse: boolean = false
): C[] {
  const n = input.length;
  if ((n & (n - 1)) !== 0) throw new Error('Input length must be a power of 2');
  const output = input.slice();
  // Bit-reversal permutation
  for (let i = 0, j = 0; i < n; i++) {
    if (i < j) [output[i], output[j]] = [output[j], output[i]];
    let m = n >> 1;
    while (m && j >= m) {
      j -= m;
      m >>= 1;
    }
    j += m;
  }
  // FFT
  for (let len = 2; len <= n; len <<= 1) {
    const wlen = field.primitiveRoot(len, inverse);
    for (let i = 0; i < n; i += len) {
      let w = field.one();
      for (let j = 0; j < len / 2; j++) {
        const u = output[i + j];
        const v = field.mul(output[i + j + len / 2], w);
        output[i + j] = field.add(u, v);
        output[i + j + len / 2] = field.sub(u, v);
        w = field.mul(w, wlen);
      }
    }
  }
  if (inverse) {
    const nInv = field.inv(field.fromNumber(n));
    for (let i = 0; i < n; i++) output[i] = field.mul(output[i], nInv);
  }
  return output;
}
