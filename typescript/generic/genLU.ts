import { IField } from './iField';

// Generic LU factorization for square matrices over any field
export function genericLU<C>(
  matrix: C[][],
  field: IField<C>
): { L: C[][]; U: C[][] } {
  const n = matrix.length;
  const L: C[][] = Array.from({ length: n }, (_, i) => Array(n).fill(field.zero()));
  const U: C[][] = Array.from({ length: n }, (_, i) => Array(n).fill(field.zero()));
  for (let i = 0; i < n; i++) {
    for (let k = i; k < n; k++) {
      let sum = field.zero();
      for (let j = 0; j < i; j++) sum = field.add(sum, field.mul(L[i][j], U[j][k]));
      U[i][k] = field.sub(matrix[i][k], sum);
    }
    for (let k = i; k < n; k++) {
      if (i === k) L[i][i] = field.one();
      else {
        let sum = field.zero();
        for (let j = 0; j < i; j++) sum = field.add(sum, field.mul(L[k][j], U[j][i]));
        L[k][i] = field.div(field.sub(matrix[k][i], sum), U[i][i]);
      }
    }
  }
  return { L, U };
}
