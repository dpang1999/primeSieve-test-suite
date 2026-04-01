use crate::helpers::lcg::Lcg;

static mut MODULUS: i32 = 2^31 -1;

fn mod_inverse(a: i32, m: i32) -> i32 {
    if a == 0 { panic! ("Inverse does not exist for zero"); }
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

fn print_matrix(a: &Vec<Vec<i32>>) {
    for row in a {
        for val in row {
            print!("{:3} ", val);
        }
        println!();
    }
    println!();
}

fn print_vector(b: &Vec<i32>) {
    for val in b {
        print!("{:3} ", val);
    }
    println!();
    println!();
}

pub fn factor(a: &mut Vec<Vec<i32>>, pivot: &mut Vec<usize>) -> i32 {
    let n = a.len();
    let m = a[0].len();
    let min_mn = std::cmp::min(m, n);
    let modulus = unsafe { MODULUS };

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
        if a[jp][j] == 0 {
            println!("Matrix is singular");
            return 1;
        }

        // Swap rows j and jp if needed
        if jp != j {
            a.swap(j, jp);
        }

        // Compute elements j+1:M of jth column
        if j < m - 1 {
            let recp = mod_inverse(a[j][j] as i32, modulus);
            for k in (j + 1)..m {
                a[k][j] = ((a[k][j] * recp)) % modulus;
            }
        }

        // Rank-1 update to trailing submatrix
        if j < min_mn - 1 {
            for ii in (j + 1)..m {
                let aii_j = a[ii][j];
                for jj in (j + 1)..n {
                    a[ii][jj] = ((a[ii][jj] - (aii_j * a[j][jj]) % modulus) + modulus) % modulus;
                }
            }
        }
    }
    0
}

pub fn solve(lu: &Vec<Vec<i32>>, pvt: &Vec<usize>, b: &mut Vec<i32>) {
    let m = lu.len();
    let n = lu[0].len();
    let modulus = unsafe { MODULUS };
    let mut ii = 0usize;
    for i in 0..m {
        let ip = pvt[i];
        let mut sum = b[ip];
        b[ip] = b[i];
        if ii == 0 {
            for j in ii..i {
            sum = ((sum - (lu[i][j] * b[j]) % modulus) + modulus) % modulus;
            }
        } else if sum == 0 {
            ii = i;
        }
        b[i] = sum;
    }

    for i in (0..n).rev() {
        let mut sum = b[i];
        for j in (i + 1)..n {
            sum = ((sum - (lu[i][j] * b[j]) % modulus) + modulus) % modulus;
        }
        b[i] = (sum * mod_inverse(lu[i][i], modulus)) % modulus;
    }
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let mut n = 4;
    if args.len() > 1 {
        n = args[1].parse().unwrap_or(4);
    }
    unsafe{MODULUS= 2_i32.pow(13)-1};
    let modulus = unsafe{MODULUS};
    
    let mut rand = Lcg::new(12345, 1345, 16645, 1013904);
    let mut a: Vec<Vec<i32>> = vec![vec![0; n]; n];
    for i in 0..n {
        let mut row_sum = 0;
        for j in 0..n {
                let val = rand.next_int() % modulus;
                a[i][j] = val;
                row_sum += val.abs();
        }
        // Set diagonal to be strictly greater than row_sum
        a[i][i] = (row_sum + rand.next_int() + 1) % modulus;
    }
   
    let mut b: Vec<i32> = (0..n).map(|_| rand.next_int() % modulus).collect();
    print_matrix(&a);
    println!("Rust specialized finite field LU");
    println!("Matrix size: {}", n);
    for i in 0..10 {
        let mut pivot: Vec<usize> = vec![0; n];
        let mut a_copy = a.clone();
        let mut b_copy = b.clone();
        factor(&mut a_copy, &mut pivot);
        solve(&a_copy, &pivot, &mut b_copy);
        println!("Iteration {} completed", i);
        /* if (i == 9) {
            println!("Final solution: ");
            print_matrix(&a_copy);
            print_vector(&b_copy);
            print_vector(&b);
        }  */
    }

}
#[allow(dead_code)]
pub fn run_algorithm() {
    main();
}
