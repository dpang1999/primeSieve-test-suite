use seeded_random::{Random,Seed};

use rust::helpers::lcg::Lcg;
use rust::helpers::prime_sieve::prime_sieve;
use rust::helpers::find_prime::{self, find_prime_congruent_one_mod_n};

pub struct FFT {}


fn mod_inverse(a: u64, m: u64) -> u64 {
    let mut a = a as i64;
    let mut m = m as i64;
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
    x1 as u64
}

fn modpow(mut base: u64, mut exp: u64, modulus: u64) -> u64 {
    if modulus == 0 { panic!("Modulus must be positive"); }
    let mut result = 1;
    base %= modulus;
    while exp > 0 {
        if exp % 2 == 1 {
            result = (result * base) % modulus;
        }
        base = (base * base) % modulus;
        exp /= 2;
    }
    result
}

fn primitive_root(modulus: u64) -> u64 {
    fn factorize(mut n: u64) -> Vec<u64> {
        let mut factors = Vec::new();
        let mut i = 2;
        while i * i <= n {
            if n % i == 0 {
                factors.push(i);
                while n % i == 0 {
                    n /= i;
                }
            }
            i += 1;
        }
        if n > 1 {
            factors.push(n);
        }
        factors
    }
    let p = modulus;
    let factors = factorize (p - 1);
    for g in 2..p {
        let mut is_root = true;
        for &factor in &factors {
            if modpow(g, (p - 1) / factor, p) == 1 {
                is_root = false;
                break;
            }
        }
        if is_root {
            return g;
        }
    }
    0

}

fn precomputeRootsOfUnity(n: usize, direction: i32, modulus: u64) -> Vec<u64> {
    let p = modulus;
    if (p - 1) % n as u64 != 0 {
        panic!("n must divide p-1 for roots of unity to exist in IntModP");
    }
    let g = primitive_root(modulus);
    println!("Primitive root: {}", g);
    let omega = modpow(g, (p-1)/n as u64, modulus);
    let mut roots = Vec::with_capacity(n as usize);
    for k in 0..n as i32 {
        let mut exponent: u32 = (k * direction % (p - 1) as i32) as u32;
        if exponent < 0 {
            exponent += (p - 1) as u32;
        }
        roots.push(modpow(omega, exponent as u64, modulus));
    }
    roots
}

impl FFT
{
    pub fn new() -> Self {
        Self {}
    }

    pub fn transform(&self, data: &mut [u64], modulus: u64, root: u64) {
        Self::transform_internal(data, -1, modulus, root);
    }

    pub fn inverse(&self, data: &mut [u64], modulus: u64, root: u64) {
        Self::transform_internal(data, 1, modulus, root);
        let nd = data.len();
        let n = nd;
        for d in 0..nd {
            data[d] = (data[d] * mod_inverse(n as u64, modulus)) % modulus;
        }
    }

    

    pub fn test(&self, data: &mut [u64], modulus: u64,root: u64) {
        let nd = data.len();
        let copy: Vec<u64> = data.iter().map(|x| *x).collect();

        self.transform(data, modulus, root);

        //println!("After transform: {}", data.iter().map(|x| x.to_string()).collect::<Vec<_>>().join(", "));
        self.inverse(data, modulus, root);
        //println!("After inverse: {}", data.iter().map(|x| x.to_string()).collect::<Vec<_>>().join(", "));
    }

    pub fn make_random(&self, n: usize) -> Vec<u64> {
        let mut rand = Lcg::new(12345, 1345, 16645, 1013904);
        let nd = 2*n;
        let mut data = Vec::with_capacity(nd);
        for _ in 0..nd {
            data.push(rand.next_int() as u64);
        }
        data
    }

  

    fn log2(n: usize) -> usize {
        let mut log = 0;
        let mut k = 1;
        while k < n {
            k *= 2;
            log += 1;
        }
        if n != (1 << log) {
            panic!("FFT: Data length is not a power of 2!: {}", n);
        }
        log
    }

    fn transform_internal(data: &mut [u64], direction: i32, modulus: u64, root: u64) {
        if data.is_empty() {
            return;
        }

        let n = data.len(); // Now n is the length of the data array (no division by 2)
        if n == 1 {
            return; // Identity operation
        }

        let logn = Self::log2(n);
        Self::bitreverse(data);

        let roots = precomputeRootsOfUnity(n, direction, modulus);
        println!("Roots: {}", roots.iter().map(|x| x.to_string()).collect::<Vec<_>>().join(", "));

        let mut dual = 1;
        for _bit in 0..logn {
            for a in 0..dual {
                let w = roots[a * (n / (2 * dual))].clone(); // Use precomputed root
                for b in (0..n).step_by(2 * dual) {
                    let i = b + a;
                    let j = b + a + dual;

                    let wd = w * &data[j] % modulus; // Twiddle factor multiplication
                    data[j] = (data[i] + modulus - wd) % modulus; // Subtract
                    data[i] = (data[i] + wd) % modulus; // Add
                }
            }
            dual *= 2;
        }
    }

    fn bitreverse(data: &mut [u64]) {
        let n = data.len();
        let nm1 = n - 1;
        let mut i = 0;
        let mut j = 0;
        while i < nm1 {
            if i < j {
                data.swap(i, j);
            }
            let mut k = n >> 1;
            while k <= j {
                j -= k;
                k >>= 1;
            }
            j += k;
            i += 1;
        }
    }

}
fn main() {
    // let mode = 0 be for testing
    let mode = 0;
    let fft = FFT::new();
    if mode != 0 {
        let args: Vec<String> = std::env::args().collect();
        let n: usize = args.get(1).and_then(|s| s.parse().ok()).unwrap_or(16);
        let modulus = find_prime_congruent_one_mod_n(n) as u64; // Example prime modulus
        let root = primitive_root(modulus);
        let mut rand = Lcg::new(12345, 1345, 16645, 1013904);
        let mut data: Vec<u64> = Vec::with_capacity(2*n);
        for _ in 0..n {
            let r = (rand.next_int() as u64) % modulus;
            data.push(r);
        }
        println!("Specialized Rust FFT Tests");
        println!("Specialized, Finite Field, n={}", n);
        for i in 0..10 {
            fft.transform(&mut data, modulus, root);
            fft.inverse(&mut data, modulus, root);
            println!("Loop {} done", i);
        }
    }
    else {
        let in1: [u64; 16] = [38, 0, 44, 87, 6, 45, 22, 93, 0, 0, 0, 0, 0, 0, 0, 0] ;
        let in2: [u64; 16] = [80, 18, 62, 90, 17, 96, 27, 97, 0, 0, 0, 0, 0, 0, 0, 0];
        //let out = [3040, 684, 5876, 11172, 5420, 16710, 12546, 20555, 16730, 15704, 21665, 5490, 13887, 4645, 9021, 0];
        let prime = 40961;
        
        let mut data1 = Vec::with_capacity(in1.len());
        let mut data2 = Vec::with_capacity(in2.len());
        for x in 0..in1.len() {
            data1.push(in1[x] as u64);
            data2.push(in2[x] as u64);
        }

        let root = primitive_root( prime );

        println!("Using modulus: {}, primitive root: {}", prime, root);
        fft.transform(&mut data1, prime, root);
        fft.transform(&mut data2, prime, root);
        
        println!("data1: {}", data1.iter().map(|x| x.to_string()).collect::<Vec<_>>().join(", "));
        println!("data2: {}", data2.iter().map(|x| x.to_string()).collect::<Vec<_>>().join(", "));

        let mut product = Vec::with_capacity(data1.len());
        for i in 0..data1.len() {
            product.push((data1[i] * data2[i]) % prime);
        }

        println!("product: {}", product.iter().map(|x| x.to_string()).collect::<Vec<_>>().join(", "));

        fft.inverse(&mut product, prime, root);
        println!("inverse product: {}", product.iter().map(|x| x.to_string()).collect::<Vec<_>>().join(", "));
       
    }
}