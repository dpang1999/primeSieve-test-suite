pub fn num_flops(m: usize, n: usize, num_iterations: usize) -> f64 {
    let md = m as f64;
    let nd = n as f64;
    let num_iter_d = num_iterations as f64;
    (md - 1.0) * (nd - 1.0) * num_iter_d * 6.0
}

pub fn execute(omega: f64, g: &mut Vec<Vec<f64>>, num_iterations: usize) {
    let m = g.len();
    let n = g[0].len();

    let omega_over_four = omega * 0.25;
    let one_minus_omega = 1.0 - omega;

    let mm1 = m - 1;
    let nm1 = n - 1;

    for _ in 0..num_iterations {
        for i in 1..mm1 {
            let gim1 = g[i - 1].clone();
            let gip1 = g[i + 1].clone();
            for j in 1..nm1 {
                g[i][j] = omega_over_four * (gim1[j] + gip1[j] + g[i][j - 1] + g[i][j + 1])
                    + one_minus_omega * g[i][j];
            }
        }
    }
}

pub fn print_matrix(a: &Vec<Vec<f64>>) {
    for row in a {
        for val in row {
            print!("{:.2} ", val);
        }
        println!();
    }
}

fn main() {
    let m = 10;
    let n = 10;
    let num_iterations = 100;
    let omega = 1.5;
    let mut g = vec![vec![0.0; n]; m];

    // Set boundary conditions
    for i in 0..m {
        g[i][0] = 0.0;         // Left edge
        g[i][n - 1] = 0.0;     // Right edge
    }
    for j in 0..n {
        g[0][j] = 100.0;       // Top edge (hot)
        g[m - 1][j] = 0.0;     // Bottom edge (cold)
    }

    println!("Initial grid:");
    print_matrix(&g);

    execute(omega, &mut g, num_iterations);

    println!("\nSteady-state temperature distribution:");
    print_matrix(&g);
}