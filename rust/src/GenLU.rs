use crate::generic::i_field::IField;
use crate::generic::double_field::DoubleField;
use crate::generic::single_field::SingleField;
use crate::generic::i_math::IMath;
use std::fmt::Display;
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
        } else if sum.coerce() == 0.0 {
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
        let mut t: U = a[j][j].abs();
        for i in (j + 1)..m {
            let ab = a[i][j].abs();
            if ab.coerce() > t.coerce() {
                jp = i;
                t = ab;
            }
        }
        pivot[j] = jp;

        // If zero pivot, factorization fails
        if a[jp][j].coerce() == 0.0 {
            return 1;
        }

        // Swap rows j and jp if needed
        if jp != j {
            a.swap(j, jp);
        }

        // Compute elements j+1:M of jth column
        if j < m - 1 {
            let recp = U::one().d(&a[j][j]);
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

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let mut n = 4;
    if args.len() > 1 {
        n = args[1].parse().unwrap_or(4);
    }

    let mut a: Vec<Vec<SingleField>> = (0..n)
        .map(|_| (0..n).map(|_| SingleField::new(rand::random::<f32>() * 1000.0)).collect())
        .collect();
    let mut b: Vec<SingleField> = (0..n)
        .map(|_| SingleField::new(rand::random::<f32>() * 1000.0))
        .collect();
    let mut pivot: Vec<usize> = vec![0; n];

    print_matrix(&a);

    factor(&mut a, &mut pivot); 

    println!("b: ");
    print_vector(&b);
    solve(&a, &pivot, &mut b);
    println!("Solution: ");
    print_vector(&b);
}