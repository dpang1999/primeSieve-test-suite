
import { IField } from './iField';
import { IMath } from './iMath';
import { LCG } from '../helpers/lcg';

// In-place LU factorization with partial pivoting, Rust-style
export function factor<C extends IField<C> & IMath<C>>(a: C[][], pivot: number[]): number {
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
      if (recp.abs() > 1000000)
      {
        console.log(`recp: ${recp.toString()}, ${recp.abs()%1_000_007}`);

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
export function solve<C extends IField<C> & IMath<C>>(lu: C[][], pvt: number[], b: C[],): void {
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
export function multiplyMatrixVector<C extends IField<C> & IMath<C>>(a: C[][], b: C[]): C[] {
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
if (require.main === module) {
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
  const rand = new LCG(12345, 1345, 65, 17);


  if (complex === 0) {
    if (mode === 1) {
      const { SingleField } = require('./singleField');
      const a = Array.from({ length: n }, () => Array.from({ length: n }, () => new SingleField(rand.nextDouble() * 1000.0)));
      const b = Array.from({ length: n }, () => new SingleField(rand.nextDouble() * 1000.0));
      const pivot = Array(n).fill(0);
      const aCopy = a.map(row => row.map(x => x.copy()));
      const bCopy = b.map(x => x.copy());
      printMatrix(a);
      printVector(b);
      factor(a, pivot);
      printMatrix(a);
      solve(a, pivot, b);
      const product = multiplyMatrixVector(aCopy, b);
      printVector(product);
    } else if (mode === 2) {
      const { DoubleField } = require('./doubleField');
      const a = Array.from({ length: n }, () => Array.from({ length: n }, () => new DoubleField(rand.nextDouble() * 1000.0)));
      const b = Array.from({ length: n }, () => new DoubleField(rand.nextDouble() * 1000.0));
      const pivot = Array(n).fill(0);
      const aCopy = a.map(row => row.map(x => x.copy()));
      const bCopy = b.map(x => x.copy());
      factor(a, pivot);
      solve(a, pivot, b);
      const product = multiplyMatrixVector(aCopy, b);
      printVector(product);
    } else {
      const { IntModP } = require('./intModP');
      const primes = primeSieve(Math.floor(rand.nextDouble() * 36340 + 10000));
      const prime = primes.findLastIndex((isComposite) => !isComposite);
      console.log(prime);
      const a = Array.from({ length: n }, () => Array.from({ length: n }, () => new IntModP(rand.nextInt(), prime)));
      const b = Array.from({ length: n }, () => new IntModP(rand.nextInt(), prime));
      const pivot = Array(n).fill(0);
      const aCopy = a.map(row => row.map(x => x.copy()));
      const bCopy = b.map(x => x.copy());
      //printMatrix(a);
      //printVector(b);
      factor(a, pivot);
      //printMatrix(a);
      solve(a, pivot, b);
      //const product = multiplyMatrixVector(aCopy, b);
      //printVector(product);
    }
  } else {
    // Complex not implemented in this stub, but could be added with a ComplexField class
    console.log('Complex field not implemented in this TypeScript stub.');
  }
}
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