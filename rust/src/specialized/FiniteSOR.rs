static mut MODULUS: i32 = 7727;

fn mod_inverse(a: i32, m: i32) -> i32 {
    if a == 0 {panic!{"Inverse does not exist"};}
    let mut m = m;
    let mut a = a;
    let (mut x0, mut x1) = (0, 1);
    let m0 = m;
    while a > 1 {
        let q = a / m;
        let t = m;
        m = a % m;
        a = t;
        let t = x0;
        x0 = x1 - q * x0;
        x1 = t;
    }
    if x1 < 0 { x1 += m0; }
    x1
}


pub fn execute(omega: i32, g: &mut Vec<Vec<i32>>, num_iterations: usize) {
    let m = g.len();
    let n = g[0].len();
    let modulus = unsafe { MODULUS };

    let omega_over_four = omega * mod_inverse(4, modulus) % modulus;
    let one_minus_omega = 1 - omega;

    let mm1 = m - 1;
    let nm1 = n - 1;

    for _ in 0..num_iterations {
        for i in 1..mm1 {
            let gim1 = g[i - 1].clone();
            let gip1 = g[i + 1].clone();
            for j in 1..nm1 {
                g[i][j] = omega_over_four * (gim1[j] + gip1[j] + g[i][j - 1] + g[i][j + 1])
                    + one_minus_omega * g[i][j];
                g[i][j] %= modulus;
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
    // arg1 = grid size n (nxn)
    let args = std::env::args().collect::<Vec<String>>();
    let n: usize = args.get(1).and_then(|s| s.parse().ok()).unwrap_or(10000);
    let m = n;
    let modulus = unsafe {MODULUS};
    let num_iterations = 10000;
    let omega = 3 * mod_inverse(2, modulus);
    let mut g = vec![vec![0; n]; m];

    // Set boundary conditions
/*     for i in 0..m {
        g[i][0] = 0.0;         // Left edge
        g[i][n - 1] = 0.0;     // Right edge
    } */
    for j in 0..n {
        g[0][j] = 100;       // Top edge (hot)
        //g[m - 1][j] = 0.0;     // Bottom edge (cold)
    }

    //println!("Initial grid:");
    //print_matrix(&g);

    println!("Rust specialized finite field SOR");
    println!("Grid size: {}x{}", m, n);
    println!("Number of iterations: {}", num_iterations);
    execute(omega, &mut g, num_iterations);

    //println!("\nSteady-state temperature distribution:");
    //print_matrix(&g);
}
#[allow(dead_code)]
pub fn run_algorithm() {
    main();
}
