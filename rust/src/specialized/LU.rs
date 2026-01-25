use rust::helpers::lcg::Lcg;

fn print_matrix(a: &Vec<Vec<f64>>) {
    for row in a {
        for val in row {
            print!("{:.3} ", val);
        }
        println!();
    }
}

fn print_vector(b: &Vec<f64>) {
    for val in b {
        print!("{:.3} ", val);
    }
    println!();
}

pub fn factor(a: &mut Vec<Vec<f64>>, pivot: &mut Vec<usize>) -> i32 {
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
        if a[jp][j] == 0.0 {
            return 1;
        }

        // Swap rows j and jp if needed
        if jp != j {
            a.swap(j, jp);
        }

        // Compute elements j+1:M of jth column
        if j < m - 1 {
            let recp = 1.0 / a[j][j];
            for k in (j + 1)..m {
                a[k][j] *= recp;
            }
        }

        // Rank-1 update to trailing submatrix
        if j < min_mn - 1 {
            for ii in (j + 1)..m {
                let aii_j = a[ii][j];
                for jj in (j + 1)..n {
                    a[ii][jj] -= aii_j * a[j][jj];
                }
            }
        }
    }
    0
}

pub fn solve(lu: &Vec<Vec<f64>>, pvt: &Vec<usize>, b: &mut Vec<f64>) {
    let m = lu.len();
    let n = lu[0].len();
    let mut ii = 0usize;
    for i in 0..m {
        let ip = pvt[i];
        let mut sum = b[ip];
        b[ip] = b[i];
        if ii == 0 {
            for j in ii..i {
                sum -= lu[i][j] * b[j];
            }
        } else if sum == 0.0 {
            ii = i;
        }
        b[i] = sum;
    }

    for i in (0..n).rev() {
        let mut sum = b[i];
        for j in (i + 1)..n {
            sum -= lu[i][j] * b[j];
        }
        b[i] = sum / lu[i][i];
    }
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let mut n = 4;
    if args.len() > 1 {
        n = args[1].parse().unwrap_or(4);
    }

    let mut rand = Lcg::new(12345, 1345, 16645, 1013904);
    let mut a: Vec<Vec<f64>> = (0..n)
        .map(|_| (0..n).map(|_| rand.next_double() * 1000.0).collect())
        .collect();
    let mut b: Vec<f64> = (0..n).map(|_| rand.next_double() * 1000.0).collect();
    let mut pivot: Vec<usize> = vec![0; n];

    //print_matrix(&a);

    factor(&mut a, &mut pivot); 

    // Only needed for debugging
    //println!("b: ");
    //print_vector(&b);
    solve(&a, &pivot, &mut b); 
    //println!("Solution (x): ");
    //print_vector(&b);
}