import { LCG } from "../helpers/lcg";

function printMatrix(a: number[][]) {
  for (const row of a) {
    console.log(row.map(val => val.toFixed(3)).join(' '));
  }
}

function printVector(b: number[]) {
  console.log(b.map(val => val.toFixed(3)).join(' '));
}

function main() {
  const n = parseInt(process.argv[2] ?? "4", 10);
  const rand = new LCG(12345, 1345, 16645, 1013904);
  const a: number[][] = Array.from({ length: n }, () =>
    Array.from({ length: n }, () => rand.nextDouble() * 1000)
  );
  const b: number[] = Array.from({ length: n }, () => rand.nextDouble() * 1000);
  const pivot: number[] = Array(n).fill(0);
  //printMatrix(a);
  factor(a, pivot);
  //console.log('b:');
  //printVector(b);
  solve(a, pivot, b);
  //console.log('Solution (x):');
  //printVector(b);
}

if (require.main === module) {
  main();
}


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
      const recp = 1.0 / a[j][j];
      for (let k = j + 1; k < m; k++) {
        a[k][j] *= recp;
      }
    }
    // Rank-1 update to trailing submatrix
    if (j < min_mn - 1) {
      for (let ii = j + 1; ii < m; ii++) {
        const aii_j = a[ii][j];
        for (let jj = j + 1; jj < n; jj++) {
          a[ii][jj] -= aii_j * a[j][jj];
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
        sum -= lu[i][j] * b[j];
      }
    } else if (sum === 0.0) {
      ii = i;
    }
    b[i] = sum;
  }
  for (let i = n - 1; i >= 0; i--) {
    let sum = b[i];
    for (let j = i + 1; j < n; j++) {
      sum -= lu[i][j] * b[j];
    }
    b[i] = sum / lu[i][i];
  }
}
