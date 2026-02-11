use crate::generic::i_field::IField;
use crate::generic::double_field::DoubleField;
use crate::generic::single_field::SingleField;
//use crate::generic::i_math::IMath;
use crate::generic::int_mod_p::IntModP;
use crate::generic::int_mod_p::set_modulus;
use crate::generic::complex_field::ComplexField;
use std::fmt::Display;
use rust::helpers::lcg::Lcg;
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
    // arg1 = grid size n (nxn)
    // arg2 = mode (1=SingleField, 2=DoubleField, else IntModP)
    // arg3 = complex_bool (0=real, 1=complex)
    let args = std::env::args().collect::<Vec<String>>();
    let n: usize = args.get(1).and_then(|s| s.parse().ok()).unwrap_or(16);
    let m = n;
    let field_type: usize = args.get(2).and_then(|s| s.parse().ok()).unwrap_or(0);
    let complex_bool: usize = args.get(3).and_then(|s| s.parse().ok()).unwrap_or(0);
    let num_iterations = 1000;
    let mut rand = Lcg::new(12345, 1345, 16645, 1013904);
    
    if complex_bool == 0 {
        println!("Not Complex");
        if field_type == 1 {
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

            //println!("Initial grid:");
            //print_matrix(&g);

            execute(omega, &mut g, num_iterations);

            //println!("\nSteady-state temperature distribution:");
            //print_matrix(&g);
        }
        else if field_type == 2 {
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

            //println!("Initial grid:");
            print_matrix(&g);

            execute(omega, &mut g, num_iterations);

            //println!("\nSteady-state temperature distribution:");
            //print_matrix(&g);
        }
        else {
            println!("Using IntModP");
            let primes = prime_sieve((rand.next_double()*36340.0+10000.0) as usize); // max i32 is 2147483647, sqrt is 46340.95 to avoid overflow
            let prime = primes.last().expect("No prime found in the range");
            set_modulus(*prime as u64);

            let omega = IntModP::new(3).d(&IntModP::new(2)); // 1.5 mod 449
            let mut g = vec![vec![omega.zero(); n]; m];

            // Set boundary conditions
            for i in 0..m {
                g[i][0] = omega.zero();         // Left edge
                g[i][n - 1] = omega.zero();     // Right edge
            }
            for j in 0..n {
                g[0][j] = IntModP::new(100);       // Top edge (hot)
                g[m - 1][j] = omega.zero();     // Bottom edge (cold)
            }

            //println!("Initial grid:");
            //print_matrix(&g);

            execute(omega, &mut g, num_iterations);

            //println!("\nSteady-state temperature distribution:");
            //print_matrix(&g);
        }
    }
    else {
        println!("Complex");
        if field_type == 1 {
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

            //println!("Initial grid:");
            //print_matrix(&g);

            execute(omega, &mut g, num_iterations);

            //println!("\nSteady-state temperature distribution:");
            //print_matrix(&g);
        }
        else if field_type == 2 {
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

            //println!("Initial grid:");
            //print_matrix(&g);

            execute(omega, &mut g, num_iterations);

            //println!("\nSteady-state temperature distribution:");
            //print_matrix(&g);
        }
        else {
            println!("Using IntModP");
            let primes = prime_sieve((rand.next_double()*36340.0+10000.0) as usize); // max i32 is 2147483647, sqrt is 46340.95 to avoid overflow
            let prime = primes.last().expect("No prime found in the range");
            set_modulus(*prime as u64);

            let omega = ComplexField::new(IntModP::new(3).d(&IntModP::new(2)),IntModP::new(0)); // 1.5 mod 449
            let mut g = vec![vec![omega.zero(); n]; m];

            // Set boundary conditions
            for i in 0..m {
                g[i][0] = omega.zero();         // Left edge
                g[i][n - 1] = omega.zero();     // Right edge
            }
            for j in 0..n {
                g[0][j] = ComplexField::new(IntModP::new(100), IntModP::new(1));       // Top edge (hot)
                g[m - 1][j] = omega.zero();     // Bottom edge (cold)
            }

            //println!("Initial grid:");
            //print_matrix(&g);

            execute(omega, &mut g, num_iterations);

            //println!("\nSteady-state temperature distribution:");
            //print_matrix(&g);
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