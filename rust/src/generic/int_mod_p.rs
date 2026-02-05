use std::fmt;
use crate::generic::i_field::IField;
use crate::generic::i_ordered::IOrdered;
use crate::generic::i_math::IMath;
use crate::generic::i_primitive_roots::IPrimitiveRoots;
use std::hash::Hash;
use std::cmp::Eq;
use std::sync::OnceLock;
#[derive(Debug)]
pub struct IntModP {
    pub i: u128,
}

pub static MODULUS: OnceLock<u128> = OnceLock::new();

fn get_modulus() -> u128 {
    *MODULUS.get().expect("Modulus not set")
}

fn mod_inverse(a: u128, p: u128) -> u128 {
    let (mut t, mut new_t) = (0 as i128, 1 as i128);
    let (mut r, mut new_r) = (p as i128, a.rem_euclid(p) as i128);
    while new_r != 0 {
        let quotient = r / new_r;
        let temp_t = t - quotient * new_t;
        t = new_t;
        new_t = temp_t;
        let temp_r = r - quotient * new_r;
        r = new_r;
        new_r = temp_r;
    }
    if r > 1 {
        panic!("No modular inverse exists for {} mod {}", a, p);
    }
    if t < 0 {
        t += p as i128;
    }
    t as u128
}

impl IntModP {
    pub fn new(i: u128) -> Self {
        let p = get_modulus();
        IntModP { i: i.rem_euclid(p) }
    }

    pub fn copy(&self) -> IntModP {
        IntModP::new(self.i)
    }

    pub fn coerce(&self, value: f64) -> IntModP {
        IntModP::new(value as u128)
    }

    pub fn coerce_to_f64(&self) -> f64 {
        self.i as f64
    }
}

impl IField for IntModP {
    fn a(&self, o: &IntModP) -> IntModP {
        IntModP::new(self.i + o.i)
    }

    fn ae(&mut self, o: &IntModP) {
        let p = get_modulus();
        self.i = (self.i + o.i).rem_euclid(p);
    }

    fn s(&self, o: &IntModP) -> IntModP {
        let p = get_modulus();
        IntModP::new(self.i + p - o.i)
    }

    fn se(&mut self, o: &IntModP) {
        let p = get_modulus();
        self.i = (self.i + p - o.i).rem_euclid(p);
    }

    fn m(&self, o: &IntModP) -> IntModP {
        IntModP::new(self.i * o.i)
    }

    fn me(&mut self, o: &IntModP) {
        let p = get_modulus();
        self.i = (self.i * o.i).rem_euclid(p);
    }

    fn d(&self, o: &IntModP) -> IntModP {
        let p = get_modulus();
        if o.i == 0 {
            panic!("Division by zero in IntModP");
        }
        else {
            let inv = mod_inverse(o.i, p);
            IntModP::new(self.i * inv)
        }
    }
    fn de(&mut self, o: &IntModP) {
        let p = get_modulus();
        if o.i == 0 {
            panic!("Division by zero in IntModP");
        }
        else {
            let inv = mod_inverse(o.i, p);
            self.i = (self.i * inv).rem_euclid(p);
        }
    }

    fn coerce_to_f64(&self) -> f64 {
        self.i as f64
    }

    fn coerce(&self, value: f64) -> IntModP {
        IntModP::new(value as u128)
    }

    fn is_zero(&self) -> bool {
        self.i == 0
    }
    fn is_one(&self) -> bool {
        self.i == 1
    }

    fn zero(&self) -> IntModP {
        IntModP::new(0)
    }
    fn one(&self) -> IntModP {
        IntModP::new(1)
    }
    fn copy(&self) -> IntModP {
        self.copy()
    }
}

impl IOrdered for IntModP {
    fn lt(&self, o: &IntModP) -> bool {
        self.i < o.i
    }

    fn le(&self, o: &IntModP) -> bool {
        self.i <= o.i
    }

    fn gt(&self, o: &IntModP) -> bool {
        self.i > o.i
    }

    fn ge(&self, o: &IntModP) -> bool {
        self.i >= o.i
    }

    fn e(&self, o: &IntModP) -> bool {
        self.i == o.i
    }
}

impl IMath for IntModP {
    fn abs(&self) -> f64 {
        self.i as f64
    }

    fn sqrt(&mut self) {
        panic!("Square root not implemented for IntModP");
    }
}

impl fmt::Display for IntModP {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "IntModP({})", self.i)
    }
}

impl Clone for IntModP {
    fn clone(&self) -> Self {
        self.copy()
    }
}

// Factorize a number
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

 fn mod_pow(base: u128, exp: u128, modulus: u128) -> u128 {
        if modulus <= 0 {
            panic!("Modulus must be positive");
        }

        let mut result = 1;
        let mut base = base % modulus;
        let mut exp = exp;

        while exp > 0 {
            if exp % 2 == 1 {
                result = (result * base) % modulus;
            }
            base = (base * base) % modulus;
            exp /= 2;
        }

        result
    }


impl IPrimitiveRoots<IntModP> for IntModP {
    fn primitive_root(&self, n: u128) -> Self {
        let p = get_modulus();
        if n == 0 || n >= p {
            panic!("n must be in range [1, p-1]");
        }
        let factors = factorize(p as u64 - 1);
        for g in 2..p {
            let mut is_root = true;
            for &factor in &factors {
                if mod_pow(g, (p - 1) / factor as u128, p) == 1 {
                    is_root = false;
                    break;
                }
            }
            if is_root {
                return Self::new(g);
            }
        }
        Self::new(0)
    }

    fn pow(&self, exp: u128) -> IntModP {
        let p = get_modulus();
        IntModP::new(mod_pow(self.i, exp as u128, p))
    }

    fn precomputeRootsOfUnity(&self, n: u32, direction: i32) -> Vec<IntModP> {
        let p = get_modulus();
        if (p - 1) % n as u128 != 0 {
            panic!("n must divide p-1 for roots of unity to exist in IntModP");
        }
        let g = self.primitive_root(p - 1);
        let omega = g.pow((p - 1) / (n as u128));
        let mut roots = Vec::with_capacity(n as usize);
        for k in 0..n as i32 {
            let mut exponent: u128 = (k * direction % (p - 1) as i32) as u128;
            if exponent < 0 {
                exponent += (p - 1) as u128;
            }
            roots.push(omega.pow(exponent));
        }
        roots
    }
}

impl Hash for IntModP {
    fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
        self.i.hash(state);
    }
}

impl PartialEq for IntModP {
    fn eq(&self, other: &Self) -> bool {
        self.i == other.i
    }
}

impl Eq for IntModP {}