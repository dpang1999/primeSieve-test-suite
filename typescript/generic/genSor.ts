
import { IField } from './iField';
import { IMath } from './iMath';
import { LCG } from '../helpers/lcg';

// SOR for grid-based problems, Rust-style
export function executeSOR<C extends IField<C>>(omega: C, g: C[][], numIterations: number): void {
  const m = g.length;
  const n = g[0].length;
  let four = g[0][0].one();
  for (let i = 0; i < 2; i++) four = four.a(four); // four = 4*one
  const omegaOverFour = omega.d(four);
  const oneMinusOmega = omega.one().s(omega);
  const mm1 = m - 1;
  const nm1 = n - 1;
  for (let iter = 0; iter < numIterations; iter++) {
    for (let i = 1; i < mm1; i++) {
      for (let j = 1; j < nm1; j++) {
        const up = g[i - 1][j];
        const down = g[i + 1][j];
        const left = g[i][j - 1];
        const right = g[i][j + 1];
        const center = g[i][j];
        const neighborSum = up.a(down).a(left.a(right));
        const newVal = omegaOverFour.m(neighborSum).a(oneMinusOmega.m(center));
        g[i][j] = newVal;
      }
    }
  }
}

// Utility: print matrix
export function printMatrix<C>(a: C[][]): void {
  for (const row of a) {
    console.log(row.map(String).join(' '));
  }
  console.log();
}
if (require.main === module) {
  // Usage: tsx genSor.ts [n] [mode] [complex]
  // mode: 1=SingleField, 2=DoubleField, else IntModP
  // complex: 0=real, 1=complex
  const args = process.argv.slice(2);
  const n = parseInt(args[0] || '16', 10);
  const m = n;
  const mode = parseInt(args[1] || '2', 10);
  const complex = parseInt(args[2] || '0', 10);
  const numIterations = 1_000_000;
  const rand = new LCG(12345, 1345, 65, 17);

  if (complex === 0) {
    if (mode === 1) {
      const { SingleField } = require('./singleField');
      const omega = new SingleField(1.5);
      const g = Array.from({ length: m }, () => Array(n).fill(omega.zero()));
      for (let i = 0; i < m; i++) {
        g[i][0] = omega.zero();
        g[i][n - 1] = omega.zero();
      }
      for (let j = 0; j < n; j++) {
        g[0][j] = new SingleField(100.0);
        g[m - 1][j] = omega.zero();
      }
      executeSOR(omega, g, numIterations);
      printMatrix(g);
    } else if (mode === 2) {
      const { DoubleField } = require('./doubleField');
      const omega = new DoubleField(1.5);
      const g = Array.from({ length: m }, () => Array(n).fill(omega.zero()));
      for (let i = 0; i < m; i++) {
        g[i][0] = omega.zero();
        g[i][n - 1] = omega.zero();
      }
      for (let j = 0; j < n; j++) {
        g[0][j] = new DoubleField(100.0);
        g[m - 1][j] = omega.zero();
      }
      executeSOR(omega, g, numIterations);
      printMatrix(g);
    } else {
      const { IntModP } = require('./intModP');
      const primes = primeSieve(Math.floor(rand.nextDouble() * 36340 + 10000))
      const prime = primes.lastIndexOf(false);
      const omega = new IntModP(3, prime).d(new IntModP(2, prime)); // 1.5 mod prime
      const g = Array.from({ length: m }, () => Array(n).fill(omega.zero()));
      for (let i = 0; i < m; i++) {
        g[i][0] = omega.zero();
        g[i][n - 1] = omega.zero();
      }
      for (let j = 0; j < n; j++) {
        g[0][j] = new IntModP(100, prime);
        g[m - 1][j] = omega.zero();
      }
      executeSOR(omega, g, numIterations);
      printMatrix(g);
    }
  } else {
    // Complex not implemented in this stub, but could be added with a ComplexField class
    console.log('Complex field not implemented in this TypeScript stub.');
  }
}
