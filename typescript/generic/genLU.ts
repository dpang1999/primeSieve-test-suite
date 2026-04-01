
import { IField } from './iField';
import { IMath } from './iMath';
import { LCG } from '../helpers/lcg.js';
import { ICopiable } from './iCopiable';
import {SingleField} from './singleField.js';
import {DoubleField} from './doubleField.js';
import {IntModP} from './intModP.js';
import { mainModule } from 'node:process';

// In-place LU factorization with partial pivoting, Rust-style
export function factor<C extends IField<C> & IMath<C> & ICopiable<C>>(a: C[][], pivot: number[]): number {
  const n = a.length;
  const m = a[0].length;
  const min_mn = Math.min(m, n);
  for (let j = 0; j < min_mn; j++) {
    // Find pivot in column j and test for singularity
    let jp = j;
    let t = a[j][j].abs();
    for (let i = j + 1; i < m; i++) {
      const ab = a[i][j].abs();
      if (ab > t) {
        jp = i;
        t = ab;
      }
    }
    pivot[j] = jp;
    if (a[jp][j].is_zero()) return 1; // singular
    if (jp !== j) {
      [a[j], a[jp]] = [a[jp], a[j]];
    }
    if (j < m - 1) {
      const recp = (a[j][j].one()).d(a[j][j]);
      if (recp.abs().coerce_to_number() > 1000000)
      {
        console.log(`recp: ${recp.toString()}, ${recp.abs().coerce_to_number()%1_000_007}`);

      }
      for (let k = j + 1; k < m; k++) {
        a[k][j] = a[k][j].m(recp);
      }
    }
    if (j < min_mn - 1) {
      for (let ii = j + 1; ii < m; ii++) {
        const aii_j = a[ii][j].copy();
        for (let jj = j + 1; jj < n; jj++) {
          a[ii][jj] = a[ii][jj].s(aii_j.m(a[j][jj]));
        }
      }
    }
  }
  return 0;
}

// Solve Ax = b using LU factorization and pivot vector
export function solve<C extends IField<C> & IMath<C> & ICopiable<C>>(lu: C[][], pvt: number[], b: C[],): void {
  const m = lu.length;
  const n = lu[0].length;
  let ii = 0;
  for (let i = 0; i < m; i++) {
    const ip = pvt[i];
    let sum = b[ip].copy();
    b[ip] = b[i].copy();
    if (ii === 0) {
      for (let j = ii; j < i; j++) {
        sum = sum.s(lu[i][j].m(b[j]));
      }
    } else if (sum.is_zero()) {
      ii = i;
    }
    b[i] = sum;
  }
  for (let i = n - 1; i >= 0; i--) {
    let sum = b[i].copy();
    for (let j = i + 1; j < n; j++) {
      sum = sum.s(lu[i][j].m(b[j]));
    }
    b[i] = sum.d(lu[i][i]);
  }
}

// Utility: multiply matrix and vector
export function multiplyMatrixVector<C extends IField<C> & IMath<C> & ICopiable<C>>(a: C[][], b: C[]): C[] {
  const m = a.length;
  const n = a[0].length;
  const product: C[] = Array(m).fill(a[0][0].zero());
  for (let i = 0; i < m; i++) {
    let sum = a[0][0].zero();
    for (let j = 0; j < n; j++) {
      sum = sum.a(a[i][j].m(b[j]));
    }
    product[i] = sum;
  }
  return product;
}

// Utility: print matrix
export function printMatrix<C>(a: C[][]): void {
  for (const row of a) {
    console.log(row.map(String).join(' '));
  }
  console.log();
}

// Utility: print vector
export function printVector<C>(b: C[]): void {
  console.log(b.map(String).join(' '));
  console.log();
}

// Main function for demonstration, similar to Rust
function main() {
  // Usage: tsx genLU.ts [n] [mode] [complex]
  // mode: 1=SingleField, 2=DoubleField, else IntModP
  // complex: 0=not complex, 1=complex
  const args = process.argv.slice(2);
  let n = 4;
  let mode = 3;
  let complex = 0;
  if (args.length > 0) n = parseInt(args[0], 10) || n;
  if (args.length > 1) mode = parseInt(args[1], 10) || mode;
  if (args.length > 2) complex = parseInt(args[2], 10) || complex;
  const rand = new LCG(12345, 1345, 16645, 1013904);


  if (complex === 0) {
    if (mode === 1) {
      console.log("TypeScript generic single field LU");
      console.log(`Matrix size: ${n}`);
      const a: SingleField[][] = Array.from({ length: n }, () =>
        Array.from({ length: n }, () => new SingleField(0))
      );
      for (let i = 0; i < n; i++) {
        let rowSum = 0;
        for (let j = 0; j < n; j++) {
          if (i !== j) {
            const val = rand.nextDouble() * 1000;
            a[i][j] = new SingleField(val);
            rowSum += Math.abs(val);
          }
        }
        // Set diagonal to be strictly greater than rowSum
        a[i][i] = new SingleField(rowSum + rand.nextDouble() * 1000 + 1);
      }
      const b = Array.from({ length: n }, () => new SingleField(rand.nextDouble() * 1000.0));
      for (let i = 0; i < 10; i++) {
        const pivot = Array(n).fill(0);
        const aCopy = a.map(row => row.map(x => x.copy()));
        const bCopy = b.map(x => x.copy());
        factor(aCopy, pivot);
        solve(aCopy, pivot, bCopy);
        console.log(`Iteration ${i} completed`);
      }
    } else if (mode === 2) {
      console.log("TypeScript generic double field LU");
      console.log(`Matrix size: ${n}`);
      const a: DoubleField[][] = Array.from({ length: n }, () =>
        Array.from({ length: n }, () => new DoubleField(0))
      );
      for (let i = 0; i < n; i++) {
        let rowSum = 0;
        for (let j = 0; j < n; j++) {
          if (i !== j) {
            const val = rand.nextDouble() * 1000;
            a[i][j] = new DoubleField(val);
            rowSum += Math.abs(val);
          }
        }
        // Set diagonal to be strictly greater than rowSum
        a[i][i] = new DoubleField(rowSum + rand.nextDouble() * 1000 + 1);
      }
      const b = Array.from({ length: n }, () => new DoubleField(rand.nextDouble() * 1000.0));
      for (let i = 0; i < 10; i++) {
        const pivot = Array(n).fill(0);
        const aCopy = a.map(row => row.map(x => x.copy()));
        const bCopy = b.map(x => x.copy());
        factor(aCopy, pivot);
        solve(aCopy, pivot, bCopy);
        console.log(`Iteration ${i} completed`);
      }
    } else {
      console.log("TypeScript generic finite field LU");
      console.log(`Matrix size: ${n}`);
      const prime = 2**13 -1
      IntModP.setModulus(prime);
      const a: IntModP[][] = Array.from({ length: n }, () =>
        Array.from({ length: n }, () => new IntModP(0))
      );
      for (let i = 0; i < n; i++) {
        let rowSum = 0;
        for (let j = 0; j < n; j++) {
          if (i !== j) {
            const val = rand.nextInt() % prime;
            a[i][j] = new IntModP(val);
            rowSum += Math.abs(val);
          }
        }
        // Set diagonal to be strictly greater than rowSum
        a[i][i] = new IntModP(rowSum + rand.nextInt() + 1);
      }
      const b = Array.from({ length: n }, () => new IntModP(rand.nextInt()));
      for (let i = 0; i < 10; i++) {
        const pivot = Array(n).fill(0);
        const aCopy = a.map(row => row.map(x => x.copy()));
        const bCopy = b.map(x => x.copy());
        factor(aCopy, pivot);
        solve(aCopy, pivot, bCopy);
        console.log(`Iteration ${i} completed`);
      }
    }
  } else {
    // Complex not implemented in this stub, but could be added with a ComplexField class
    console.log('Complex field not implemented in this TypeScript stub.');
  }
}
main();
function primeSieve(num:number): boolean[] {
    let primes = new Array(num)
    primes[0] = true;
    primes[1] = true;
    for(let i = 2; i<=num; i++) {
        if(!primes[i]) {
            let j = i;
            while (i * j < num) {
                primes[i*j] = true;
                j++;
            }
        }
    }
    return primes;
}