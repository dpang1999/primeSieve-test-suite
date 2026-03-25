import { LCG } from "../helpers/lcg.js";

let modulus = 2**19 - 1; // Mersenne prime for 32-bit integers
function modInverse(a: number, m: number): number {
  let m0 = m, x0 = 0, x1 = 1;
  if (m === 1) return 0;
  a = ((a % m) + m) % m;
  while (a > 1) {
    const q = Math.floor(a / m);
    [a, m] = [m, a % m];
    [x0, x1] = [x1 - q * x0, x0];
  }
  if (x1 < 0) x1 += m0;
  return x1;
}

function printMatrix(a: number[][]) {
  for (const row of a) {
    console.log(row.map(val => val.toFixed(3)).join(' '));
  }
}

function printVector(b: number[]) {
  console.log(b.map(val => val.toFixed(3)).join(' '));
}

function main() {
  console.log("TypeScript specialized finite field LU")
  const n = parseInt(process.argv[2] ?? "4", 10);
  console.log(`Matrix size: ${n}`);
  const rand = new LCG(987654321, 2**31 - 1, 16645, 1013904);
  const a: number[][] = Array.from({ length: n }, () => Array(n).fill(0));
  for (let i = 0; i < n; i++) {
    let rowSum = 0;
    for (let j = 0; j < n; j++) {
      if (i !== j) {
        const val = rand.nextInt();
        a[i][j] = val;
        rowSum += Math.abs(val);
      }
    }
    // Set diagonal to be strictly greater than rowSum
    a[i][i] = rowSum + rand.nextInt() + 1;
  }
  const b: number[] = Array.from({ length: n }, () => rand.nextInt());
  for (let i = 0; i < 10; i++) {
    const pivot: number[] = Array(n).fill(0);
    const aCopy = a.map(row => row.slice());
    const bCopy = b.slice();
    factor(aCopy, pivot);
    solve(aCopy, pivot, bCopy);
    console.log(`Iteration ${i} completed`);
  }
}

main();


// In-place LU factorization with pivoting, matching Rust style
export function factor(a: number[][], pivot: number[]): number {
  const n = a.length;
  const m = a[0].length;
  const min_mn = Math.min(m, n);
  for (let j = 0; j < min_mn; j++) {
    // Find pivot in column j and test for singularity
    let jp = j;
    let t = Math.abs(a[j][j]);
    for (let i = j + 1; i < m; i++) {
      const ab = Math.abs(a[i][j]);
      if (ab > t) {
        jp = i;
        t = ab;
      }
    }
    pivot[j] = jp;
    if (a[jp][j] === 0) {
      return 1;
    }
    // Swap rows j and jp if needed
    if (jp !== j) {
      const tmp = a[j];
      a[j] = a[jp];
      a[jp] = tmp;
    }
    // Compute elements j+1:M of jth column
    if (j < m - 1) {
      const recp = modInverse(a[j][j], modulus);
      for (let k = j + 1; k < m; k++) {
        a[k][j] = (a[k][j] * recp + modulus) % modulus;
      }
    }
    // Rank-1 update to trailing submatrix
    if (j < min_mn - 1) {
      for (let ii = j + 1; ii < m; ii++) {
        const aii_j = a[ii][j];
        for (let jj = j + 1; jj < n; jj++) {
          a[ii][jj] = ((a[ii][jj] - (aii_j * a[j][jj]) % modulus) + modulus) % modulus;
        }
      }
    }
  }
  return 0;
}

// Solve Ax = b using LU and pivot
export function solve(lu: number[][], pvt: number[], b: number[]): void {
  const m = lu.length;
  const n = lu[0].length;
  let ii = 0;
  for (let i = 0; i < m; i++) {
    const ip = pvt[i];
    let sum = b[ip];
    b[ip] = b[i];
    if (ii === 0) {
      for (let j = ii; j < i; j++) {
        sum = (sum - (lu[i][j] * b[j]) % modulus + modulus) % modulus;
      }
    } else if (sum === 0.0) {
      ii = i;
    }
    b[i] = sum;
  }
  for (let i = n - 1; i >= 0; i--) {
    let sum = b[i];
    for (let j = i + 1; j < n; j++) {
      sum = (sum - (lu[i][j] * b[j]) % modulus + modulus) % modulus;
    }
    b[i] = ((sum * modInverse(lu[i][i], modulus)) + modulus) % modulus;
  }
}
