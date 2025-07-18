use std::fmt;
use crate::generic::i_field::IField;
use crate::generic::i_ordered::IOrdered;
use crate::generic::i_math::IMath;

pub struct IntModP {
    pub i: i32,
    pub p: i32,
    pub print_short: bool,
}

fn mod_inverse(a: i32, p: i32) -> i32 {
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
    pub fn new(i: i32, p: i32) -> Self {
        IntModP { i: i.rem_euclid(p), p, print_short: true }
    }

    pub fn copy(&self) -> IntModP {
        IntModP::new(self.i, self.p)
    }

    pub fn coerce(&self) -> i32 {
        self.i
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
 
    fn coerce(&self, value: f64) -> IntModP {
        IntModP::new(value as i32, self.p)
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