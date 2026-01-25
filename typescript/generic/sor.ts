import { IField } from './iField';

// Generic Successive Over-Relaxation (SOR) for linear systems Ax = b
type SORResult<C> = { x: C[]; iterations: number };

export function genericSOR<C>(
  A: C[][],
  b: C[],
  field: IField<C>,
  omega: number = 1.0,
  tol: number = 1e-8,
  maxIter: number = 1000
): SORResult<C> {
  const n = b.length;
  let x = Array(n).fill(field.zero());
  for (let iter = 0; iter < maxIter; iter++) {
    let maxDiff = 0;
    for (let i = 0; i < n; i++) {
      let sigma = field.zero();
      for (let j = 0; j < n; j++) {
        if (j !== i) sigma = field.add(sigma, field.mul(A[i][j], x[j]));
      }
      const old = x[i];
      const aii = A[i][i];
      const bi = b[i];
      const newXi = field.add(
        field.mul(field.fromNumber(omega), field.div(field.sub(bi, sigma), aii)),
        field.mul(field.fromNumber(1 - omega), old)
      );
      x[i] = newXi;
      // @ts-ignore
      const diff = Math.abs((field.sub(newXi, old) as any).coerceToFloat?.() ?? 0);
      if (diff > maxDiff) maxDiff = diff;
    }
    if (maxDiff < tol) return { x, iterations: iter + 1 };
  }
  return { x, iterations: maxIter };
}
