use crate::generic::i_field::IField;
use crate::generic::double_field::DoubleField;
use crate::generic::single_field::SingleField;
//use crate::generic::i_math::IMath;
use crate::generic::int_mod_p::IntModP;
use crate::generic::complex_field::ComplexField;
use std::fmt::Display;
use rand::Rng;
pub mod generic;

pub fn num_flops(m: usize, n: usize, num_iterations: usize) -> f64 {
    let md = m as f64;
    let nd = n as f64;
    let num_iter_d = num_iterations as f64;
    (md - 1.0) * (nd - 1.0) * num_iter_d * 6.0
}

pub fn execute<U: IField + Display>(omega: U, g: &mut Vec<Vec<U>>, num_iterations: usize) {
    let m = g.len();
    let n = g[0].len();

    let mut four = U::one(&omega);
    for i in 0..2 {
        four = four.a(&four); // The dumbest way to make four
    }
    let omega_over_four = omega.d(&four);
    let one_minus_omega = U::one(&omega).s(&omega);

    let mm1 = m - 1;
    let nm1 = n - 1;

    for _ in 0..num_iterations {
        for i in 1..mm1 {
            for j in 1..nm1 {
            let up    = g[i - 1][j].copy();
            let down  = g[i + 1][j].copy();
            let left  = g[i][j - 1].copy();
            let right = g[i][j + 1].copy();
            let center = g[i][j].copy();

            let neighbor_sum = up.a(&down).a(&left).a(&right);
            let new_val = omega_over_four.m(&neighbor_sum)
                .a(&one_minus_omega.m(&center));
            g[i][j] = new_val;
            }
        }
    }
}

fn print_matrix<T: Display>(a: &Vec<Vec<T>>) {
    for row in a {
        for val in row {
            print!("{} ", val);
        }
        println!();
    }
    println!();
}

fn main() {
    let args = std::env::args().collect::<Vec<String>>();
    let mut m = 10;
    let mut n = 10;
    let mut mode = 2;
    let mut complex_bool = 0;
    let num_iterations = 100;
    let mut rng = rand::rng();
    if args.len() > 1 {
        n = args[1].parse().unwrap_or(4);
        m = n;
    }
    if args.len() > 2 {
        mode = args[2].parse().unwrap_or(4);
    }
    if args.len() > 3 {
        complex_bool = args[3].parse().unwrap_or(0);
    }
    
    if complex_bool == 0 {
        println!("Not Complex");
        if mode == 1 {
            println!("Using SingleField");
            let omega = SingleField::new(1.5);
            let mut g = vec![vec![omega.zero(); n]; m];

            // Set boundary conditions
            for i in 0..m {
                g[i][0] = omega.zero();         // Left edge
                g[i][n - 1] = omega.zero();     // Right edge
            }
            for j in 0..n {
                g[0][j] = SingleField::new(100.0);       // Top edge (hot)
                g[m - 1][j] = omega.zero();     // Bottom edge (cold)
            }

            println!("Initial grid:");
            print_matrix(&g);

            execute(omega, &mut g, num_iterations);

            println!("\nSteady-state temperature distribution:");
            print_matrix(&g);
        }
        else if mode == 2 {
            println!("Using DoubleField");
            let omega = DoubleField::new(1.5);
            let mut g = vec![vec![omega.zero(); n]; m];

            // Set boundary conditions
            for i in 0..m {
                g[i][0] = omega.zero();         // Left edge
                g[i][n - 1] = omega.zero();     // Right edge
            }
            for j in 0..n {
                g[0][j] = DoubleField::new(100.0);       // Top edge (hot)
                g[m - 1][j] = omega.zero();     // Bottom edge (cold)
            }

            println!("Initial grid:");
            print_matrix(&g);

            execute(omega, &mut g, num_iterations);

            println!("\nSteady-state temperature distribution:");
            print_matrix(&g);
        }
        else {
            println!("Using IntModP");
            let primes = prime_sieve(rng.random_range(10000..46340)); // max i32 is 2147483647, sqrt is 46340.95 to avoid overflow
            let prime = primes.last().expect("No prime found in the range");

            let omega = IntModP::new(3, *prime as u128).d(&IntModP::new(2, *prime as u128)); // 1.5 mod 449
            let mut g = vec![vec![omega.zero(); n]; m];

            // Set boundary conditions
            for i in 0..m {
                g[i][0] = omega.zero();         // Left edge
                g[i][n - 1] = omega.zero();     // Right edge
            }
            for j in 0..n {
                g[0][j] = IntModP::new(100, *prime as u128);       // Top edge (hot)
                g[m - 1][j] = omega.zero();     // Bottom edge (cold)
            }

            println!("Initial grid:");
            print_matrix(&g);

            execute(omega, &mut g, num_iterations);

            println!("\nSteady-state temperature distribution:");
            print_matrix(&g);
        }
    }
    else {
        println!("Complex");
        if mode == 1 {
            println!("Using SingleField");
            let omega = ComplexField::new(SingleField::new(1.5), SingleField::new(0.0));
            let mut g = vec![vec![omega.zero(); n]; m];

            // Set boundary conditions
            for i in 0..m {
                g[i][0] = omega.zero();         // Left edge
                g[i][n - 1] = omega.zero();     // Right edge
            }
            for j in 0..n {
                g[0][j] = ComplexField::new(SingleField::new(100.0), SingleField::new(1.0));       // Top edge (hot)
                g[m - 1][j] = omega.zero();     // Bottom edge (cold)
            }

            println!("Initial grid:");
            print_matrix(&g);

            execute(omega, &mut g, num_iterations);

            println!("\nSteady-state temperature distribution:");
            print_matrix(&g);
        }
        else if mode == 2 {
            println!("Using DoubleField");
            let omega = ComplexField::new(DoubleField::new(1.5), DoubleField::new(0.0));
            let mut g = vec![vec![omega.zero(); n]; m];

            // Set boundary conditions
            for i in 0..m {
                g[i][0] = omega.zero();         // Left edge
                g[i][n - 1] = omega.zero();     // Right edge
            }
            for j in 0..n {
                g[0][j] = ComplexField::new(DoubleField::new(100.0), DoubleField::new(1.0));       // Top edge (hot)
                g[m - 1][j] = omega.zero();     // Bottom edge (cold)
            }

            println!("Initial grid:");
            print_matrix(&g);

            execute(omega, &mut g, num_iterations);

            println!("\nSteady-state temperature distribution:");
            print_matrix(&g);
        }
        else {
            println!("Using IntModP");
            let primes = prime_sieve(rng.random_range(10000..46340)); // max i32 is 2147483647, sqrt is 46340.95 to avoid overflow
            let prime = primes.last().expect("No prime found in the range");

            let omega = ComplexField::new(IntModP::new(3, *prime as u128).d(&IntModP::new(2, *prime as u128)),IntModP::new(0,*prime as u128)); // 1.5 mod 449
            let mut g = vec![vec![omega.zero(); n]; m];

            // Set boundary conditions
            for i in 0..m {
                g[i][0] = omega.zero();         // Left edge
                g[i][n - 1] = omega.zero();     // Right edge
            }
            for j in 0..n {
                g[0][j] = ComplexField::new(IntModP::new(100, *prime as u128), IntModP::new(1,*prime as u128));       // Top edge (hot)
                g[m - 1][j] = omega.zero();     // Bottom edge (cold)
            }

            println!("Initial grid:");
            print_matrix(&g);

            execute(omega, &mut g, num_iterations);

            println!("\nSteady-state temperature distribution:");
            print_matrix(&g);
        }
    }
}

fn prime_sieve(num:usize) -> Vec<i32> {
    let mut numbers:Vec<bool> = vec![true;num];
    let mut prime_numbers:Vec<i32> = vec![];
    numbers[0] = false;
    numbers[1] = false;
    for i in 2..num {
        //interestingly i never takes the value of num, non-inclusive end range
        if numbers[i] {
            prime_numbers.push(i as i32);
            let mut j: usize = i;
            let mut current: usize = j*i;
            while current<num {
                numbers[current] = false;
                j+=1;
                current = j*i;
            }
        
        }
    }
    prime_numbers
}