use crate::generic::i_field::IField;
use crate::generic::double_field::DoubleField;
use crate::generic::single_field::SingleField;
use crate::generic::i_math::IMath;
use crate::generic::int_mod_p::IntModP;
use crate::generic::int_mod_p::MODULUS;
use crate::generic::complex_field::ComplexField;
use std::fmt::Display;
use rust::helpers::lcg::Lcg;
pub mod generic;

pub fn solve<U: IField>(
    lu: &Vec<Vec<U>>,
    pvt: &Vec<usize>,
    b: &mut Vec<U>,
) {
    let m = lu.len();
    let n = lu[0].len();
    let mut ii = 0usize;

    for i in 0..m {
        let ip = pvt[i];
        let mut sum = b[ip].copy();

        b[ip] = b[i].copy();
        if ii == 0 {
            for j in ii..i {
                sum = sum.s(&lu[i][j].m(&b[j]));
            }
        } else if sum.is_zero() {
            ii = i;
        }
        b[i] = sum;
    }

    for i in (0..n).rev() {
        let mut sum = b[i].copy();
        for j in (i + 1)..n {
            sum = sum.s(&lu[i][j].m(&b[j]));
        }
        b[i] = sum.d(&lu[i][i]);
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

fn print_vector<T: Display>(b: &Vec<T>) {
    for val in b {
        print!("{} ", val);
    }
    println!();
    println!();
}



pub fn factor<U: IField + IMath>(
    a: &mut Vec<Vec<U>>,
    pivot: &mut Vec<usize>,
) -> i32 {
    let n = a.len();
    let m = a[0].len();
    let min_mn = std::cmp::min(m, n);

    for j in 0..min_mn {
        // Find pivot in column j and test for singularity
        let mut jp = j;
        let mut t = a[j][j].abs();
        for i in (j + 1)..m {
            let ab = a[i][j].abs();
            if ab > t {
                jp = i;
                t = ab;
            }
        }
        pivot[j] = jp;

        // If zero pivot, factorization fails
        if a[jp][j].is_zero() {
            return 1;
        }

        // Swap rows j and jp if needed
        if jp != j {
            a.swap(j, jp);
        }

        // Compute elements j+1:M of jth column
        if j < m - 1 {
            let recp = a[j][j].one().d(&a[j][j]);
            for k in (j + 1)..m {
                a[k][j] = a[k][j].m(&recp);
            }
        }

        // Rank-1 update to trailing submatrix
        if j < min_mn - 1 {
            for ii in (j + 1)..m {
                let aii_j = a[ii][j].copy();
                for jj in (j + 1)..n {
                    a[ii][jj] = a[ii][jj].s(&aii_j.m(&a[j][jj]));
                }
            }
        }
    }
    0
}

fn run<T: IField + IMath + Clone + Display>(
    mut a: Vec<Vec<T>>,
    mut b: Vec<T>,
    mut pivot: Vec<usize>,
) {
    print_matrix(&a);
    let a_copy = a.clone();
    factor(&mut a, &mut pivot);
    //print_matrix(&a);
    println!("b: ");
    print_vector(&b);
    let b_copy = b.clone();
    solve(&a, &pivot, &mut b);
    println!("Solution: ");
    print_vector(&b);
    let product = multiplyMatrices(a_copy, b);
    //print_vector(&product);

    // RMS diff between b_copy and product
    /*let mut rms_diff = 0.0;
    for i in 0..b_copy.len() {
        let diff = b_copy[i].s(&product[i]);
        let diff_f64 = diff.coerce_to_f64();
        rms_diff += diff_f64 * diff_f64;
    }
    rms_diff = (rms_diff / (b_copy.len() as f64)).sqrt();
    println!("RMS difference between original b and A*x: {}", rms_diff);*/
   
}

fn main() {
    // arg1 = n (matrix size)
    // arg2 = field (1=SingleField, 2=DoubleField, else=int mod p)
    // arg3 = complex_bool (0=not complex, 1=complex)
    let args: Vec<String> = std::env::args().collect();
    
    let mut n = 4;
    let mut field = 3;
    let mut rand = Lcg::new(12345,1345,65,17);
    let mut complex_bool = 0;
    if args.len() > 1 {
        n = args[1].parse().unwrap_or(4);
    }
    if args.len() > 2 {
        field = args[2].parse().unwrap_or(1);
    }
    if args.len() > 3 {
        complex_bool = args[3].parse().unwrap_or(0);
    }
    
    if complex_bool == 0 {
        println!("Not Complex");
        if field == 1 {
            println!("Using SingleField");
            let a: Vec<Vec<SingleField>> = (0..n)
                .map(|_| (0..n).map(|_| SingleField::new((rand.next_double() * 1000.0) as f32)).collect())
                .collect();
            let b: Vec<SingleField> = (0..n)
                .map(|_| SingleField::new((rand.next_double() * 1000.0) as f32))
                .collect();
            let pivot: Vec<usize> = vec![0; n];
            run(a, b, pivot);
        } else if field == 2 {
            println!("Using DoubleField");
            let a: Vec<Vec<DoubleField>> = (0..n)
                .map(|_| (0..n).map(|_| DoubleField::new(rand.next_double() * 1000.0)).collect())
                .collect();
            let b: Vec<DoubleField> = (0..n)
                .map(|_| DoubleField::new(rand.next_double() * 1000.0))
                .collect();
            let pivot: Vec<usize> = vec![0; n];
            run(a, b, pivot);
        } else {
            println!("Using int mod p");
            let primes = prime_sieve((rand.next_double()*36340.0+10000.0) as usize); // max i32 is 2147483647, sqrt is 46340.95 to avoid overflow
            let prime = primes.last()
                .expect("No prime found in the range");
            MODULUS.set(*prime as u128).unwrap();
            let a: Vec<Vec<IntModP>> = (0..n)
                .map(|_| (0..n).map(|_| IntModP::new(rand.next_int() as u128)).collect())
                .collect();
            let b: Vec<IntModP> = (0..n)
                .map(|_| IntModP::new(rand.next_int() as u128))
                .collect();
            let pivot: Vec<usize> = vec![0; n];
            run(a, b, pivot);
        }
    }
    else if (complex_bool == 1){
        println!("Complex");
        if field == 1 {
            println!("Using SingleField");
            let a: Vec<Vec<SingleField>> = (0..n)
                .map(|_| (0..n).map(|_| SingleField::new((rand.next_double() * 1000.0) as f32)).collect())
                .collect();
            let b: Vec<SingleField> = (0..n)
                .map(|_| SingleField::new((rand.next_double() * 1000.0) as f32))
                .collect();
            let pivot: Vec<usize> = vec![0; n];
            run(a, b, pivot);
        } else if field == 2 {
            println!("Using DoubleField");
            let a: Vec<Vec<ComplexField<DoubleField>>> = (0..n)
                .map(|_| (0..n)
                    .map(|_| ComplexField::new(
                        DoubleField::new(rand.next_double() * 1000.0),
                        DoubleField::new(rand.next_double() * 1000.0)
                    ))
                    .collect()
                )
                .collect(); 
            let b: Vec<ComplexField<DoubleField>> = (0..n)
                .map(|_| ComplexField::new(
                    DoubleField::new(rand.next_double() * 1000.0),
                    DoubleField::new(rand.next_double() * 1000.0)
                ))
                .collect(); 
            let pivot: Vec<usize> = vec![0; n];
            run(a, b, pivot);
        } else {
            println!("Using int mod p");
            let primes = prime_sieve((rand.next_double()*36340.0+10000.0) as usize); // max i32 is 2147483647, sqrt is 46340.95 to avoid overflow
            let prime = primes.last()
                .expect("No prime found in the range");
            MODULUS.set(*prime as u128).unwrap();
            let a: Vec<Vec<IntModP>> = (0..n)
                .map(|_| (0..n).map(|_| IntModP::new(rand.next_int() as u128)).collect())
                .collect();
            let b: Vec<IntModP> = (0..n)
                .map(|_| IntModP::new(rand.next_int() as u128))
                .collect();
            let pivot: Vec<usize> = vec![0; n];
            run(a, b, pivot);
        }
    } else {
          let a: Vec<Vec<DoubleField>> = vec![
                vec![
                    DoubleField::new(1.0),
                    DoubleField::new(2.0),
                    DoubleField::new(3.0),
                    DoubleField::new(4.0),
                ],
                vec![
                    DoubleField::new(2.0),
                    DoubleField::new(3.0),
                    DoubleField::new(4.0),
                    DoubleField::new(5.0),
                ],
                vec![
                    DoubleField::new(3.0),
                    DoubleField::new(4.0),
                    DoubleField::new(5.0),
                    DoubleField::new(6.0),
                ],
                vec![
                    DoubleField::new(4.0),
                    DoubleField::new(5.0),
                    DoubleField::new(6.0),
                    DoubleField::new(7.0),
                ],
            ];
            let b: Vec<DoubleField> = vec![
                DoubleField::new(1.0),
                DoubleField::new(2.0),
                DoubleField::new(3.0),
                DoubleField::new(4.0),
            ];
            let pivot: Vec<usize> = vec![0; 4];
            run(a, b, pivot);
        println!("Invalid complex_bool value. Use 0 for non-complex and 1 for complex.");
    }
    


}

fn multiplyMatrices<U: IField + IMath + Display + Clone>(a: Vec<Vec<U>>, b: Vec<U>) -> Vec<U> {
    let m = a.len();
    let n = a[0].len();
    let mut product: Vec<U> = vec![a[0][0].zero(); m];

    for i in 0..m {
        let mut sum = a[0][0].zero();
        for j in 0..n {
            sum = sum.a(&a[i][j].m(&b[j]));
        }
        product[i] = sum;
    }

    product
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