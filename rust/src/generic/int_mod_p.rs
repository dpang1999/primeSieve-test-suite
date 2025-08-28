use std::fmt;
use crate::generic::i_field::IField;
use crate::generic::i_ordered::IOrdered;
use crate::generic::i_math::IMath;
use crate::generic::i_primitive_roots::IPrimitiveRoots;

pub struct IntModP {
    pub i: u128,
    pub p: u128,
    pub print_short: bool,
}

fn mod_inverse(a: u128, p: u128) -> u128 {
    let (mut t, mut new_t) = (0, 1);
    let (mut r, mut new_r) = (p, a.rem_euclid(p));
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
        t += p;
    }
    t
}

impl IntModP {
    pub fn new(i: u128, p: u128) -> Self {
        IntModP { i: i.rem_euclid(p), p, print_short: true }
    }

    pub fn copy(&self) -> IntModP {
        IntModP::new(self.i, self.p)
    }

    pub fn coerce(&self, value: f64) -> IntModP {
        IntModP::new(value as u128, self.p)
    }

    pub fn coerce_to_f64(&self) -> f64 {
        self.i as f64
    }


     
}

impl IField for IntModP {
    fn a(&self, o: &IntModP) -> IntModP {
        IntModP::new(self.i + o.i, self.p)
    }

    fn ae(&mut self, o: &IntModP) {
        self.i = (self.i + o.i).rem_euclid(self.p);
    }

    fn s(&self, o: &IntModP) -> IntModP {
        IntModP::new(self.i - o.i + self.p, self.p)
    }

    fn se(&mut self, o: &IntModP) {
        self.i = (self.i - o.i + self.p).rem_euclid(self.p);
    }

    fn m(&self, o: &IntModP) -> IntModP {
        IntModP::new(self.i * o.i, self.p)
    }

    fn me(&mut self, o: &IntModP) {
        self.i = (self.i * o.i).rem_euclid(self.p);
    }

    fn d(&self, o: &IntModP) -> IntModP {
        if o.i == 0 {
            panic!("Division by zero in IntModP");
        }
        else {
            let inv = mod_inverse(o.i, self.p);
            IntModP::new(self.i * inv, self.p)
        }
       
    }
    fn de(&mut self, o: &IntModP) {
        if o.i == 0 {
            panic!("Division by zero in IntModP");
        }
        else {
            let inv = mod_inverse(o.i, self.p);
            self.i = (self.i * inv).rem_euclid(self.p);
        }
        
    }

    fn coerce_to_f64(&self) -> f64 {
        self.i as f64
    }
 
    fn coerce(&self, value: f64) -> IntModP {
        IntModP::new(value as u128, self.p)
    }

    fn is_zero(&self) -> bool {
        self.i == 0
    }
    fn is_one(&self) -> bool {
        self.i == 1
    }

    fn zero(&self) -> IntModP {
        IntModP::new(0, self.p)
    }
    fn one(&self) -> IntModP {
        IntModP::new(1, self.p)
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

    fn eq(&self, o: &IntModP) -> bool {
        self.i == o.i && self.p == o.p
    }


}

impl IMath for IntModP {
    fn abs(&self) -> f64 {
        self.i as f64
    }

    fn sqrt(&mut self) {
        // Square root in modular arithmetic is not straightforward and is not implemented here.
        panic!("Square root not implemented for IntModP");
    }
}

impl fmt::Display for IntModP {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if self.print_short {
            write!(f, "{} (mod {})", self.i, self.p)
        } else {
            write!(f, "IntModP({}, {})", self.i, self.p)
        }
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
    // Primitive Root
    fn primitive_root(&self, n: u64) -> Self {
        if n == 0 || n >= self.p as u64 {
            panic!("n must be in range [1, p-1]");
        }

        // Iterate through potential primitive roots
        for g in 2..self.p {
            let mut is_root = true;

            // Check if g^k mod p != 1 for 0 < k < n
            for k in 1..n {
                if mod_pow(g, k as u128, self.p) == 1 {
                    is_root = false;
                    break;
                }
            }

            // If g is a primitive root, return it
            if is_root {
                return Self::new(g, self.p);
            }
        }

        Self::new(0, self.p) // No primitive root found
    }

   

    

    fn pow(&self, exp: i32) -> IntModP {
        IntModP::new(mod_pow(self.i, exp as u128, self.p), self.p)
    }

    fn precomputeRootsOfUnity(&self, n: u64, direction: i32) -> Vec<IntModP> {
        // This method is not implemented for IntModP
        panic!("precomputeRootsOfUnity not implemented for IntModP");
    }
}