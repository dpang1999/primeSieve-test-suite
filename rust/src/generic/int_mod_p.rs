use std::fmt;
use crate::generic::i_field::IField;
use crate::generic::i_ordered::IOrdered;

pub struct IntModP {
    pub i: i32,
    pub p: i32,
    pub print_short: bool,
}

impl IntModP {
    pub fn new(i: i32, p: i32) -> Self {
        IntModP { i, p, print_short: true }
    }

    pub fn copy(&self) -> IntModP {
        IntModP::new(self.i, self.p)
    }
}

impl IField for IntModP {
    fn a(&self, o: &IntModP) -> IntModP {
        IntModP::new((self.i + o.i) % self.p, self.p)
    }

    fn ae(&mut self, o: &IntModP) {
        self.i = (self.i + o.i) % self.p;
    }

    fn s(&self, o: &IntModP) -> IntModP {
        IntModP::new((self.i - o.i + self.p) % self.p, self.p)
    }

    fn se(&mut self, o: &IntModP) {
        self.i = (self.i - o.i + self.p) % self.p;
    }

    fn m(&self, o: &IntModP) -> IntModP {
        IntModP::new((self.i * o.i) % self.p, self.p)
    }

    fn me(&mut self, o: &IntModP) {
        self.i = (self.i * o.i) % self.p;
    }

    //division in int mod p is not rigorous so I don't want to do it rn
    fn d(&self, o: &IntModP) -> IntModP {
        if o.i == 0 {
            panic!("Division by zero in IntModP");
        }
        //let inv = mod_inverse(o.i, self.p);
        //IntModP::new((self.i * inv) % self.p, self.p)
    }
    fn de(&mut self, o: &IntModP) {
        if o.i == 0 {
            panic!("Division by zero in IntModP");
        }
        //let inv = mod_inverse(o.i, self.p);
        //self.i = (self.i * inv) % self.p;
    }

    fn coerce(&self) -> f64 {
        self.i as f64
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