

// Specialized Successive Over-Relaxation (SOR) for number[][]

// In-place SOR for heat distribution (2D grid), matching Rust style
export function sorHeat(
  g: number[][],
  omega: number,
  numIterations: number
): void {
  const m = g.length;
  const n = g[0].length;
  const omega_over_four = omega * 0.25;
  const one_minus_omega = 1.0 - omega;
  const mm1 = m - 1;
  const nm1 = n - 1;
  for (let iter = 0; iter < numIterations; iter++) {
    for (let i = 1; i < mm1; i++) {
      const gim1 = g[i - 1];
      const gip1 = g[i + 1];
      for (let j = 1; j < nm1; j++) {
        g[i][j] = omega_over_four * (gim1[j] + gip1[j] + g[i][j - 1] + g[i][j + 1])
          + one_minus_omega * g[i][j];
      }
    }
  }
}


function main() {
  // arg 1 = grid size n (nxn)
  const n = parseInt(process.argv[2]) || 16;
  const m = n;
  const num_iterations = 1000;
  const omega = 1.5;
  // make m x n grid, all 0
  const g: number[][] = Array.from({ length: m }, () => Array(n).fill(0));
  // Set boundary conditions
  for (let i = 0; i < m; i++) {
    g[i][0] = 0.0;         // Left edge
    g[i][n - 1] = 0.0;     // Right edge
  }
  for (let j = 0; j < n; j++) {
    g[0][j] = 100.0;       // Top edge (hot)
    g[m - 1][j] = 0.0;     // Bottom edge (cold)
  }
  sorHeat(g, omega, num_iterations);
  // Optionally print the grid
  //for (const row of g) console.log(row.map(v => v.toFixed(2)).join(' '));
}

if (require.main === module) {
  main();
}

