use std::fmt;
use crate::generic::i_field::IField;
use crate::generic::i_math::IMath;
use crate::generic::i_ordered::IOrdered;
use crate::generic::i_trigonometric::ITrigonometric;
use std::hash::Hash;
use std::cmp::Eq;
#[derive(Debug)]
    pub struct DoubleField {
    pub d: f64,
    pub print_short: bool,
}

impl DoubleField {
    //pub static mut FCOUNT: i32 = 0;

    pub fn new(d: f64) -> Self {
        DoubleField { d, print_short: false }
    }

    pub fn copy(&self) -> DoubleField {
        DoubleField::new(self.d)
    }

    /* pub fn new_array(size: usize) -> Vec<DoubleField> {
        vec![DoubleField::new(0.0); size]
    } */
}


impl IField for DoubleField {
    fn a(&self, o: &DoubleField) -> DoubleField {
        //unsafe { DoubleField::FCOUNT += 1; }
        DoubleField::new(self.d + o.d)
    }

    fn ae(&mut self, o: &DoubleField) {
        //unsafe { DoubleField::FCOUNT += 1; }
        self.d += o.d;

    }

    fn s(&self, o: &DoubleField) -> DoubleField {
        //unsafe { DoubleField::FCOUNT += 1; }
        DoubleField::new(self.d - o.d)
    }

    fn se(&mut self, o: &DoubleField) {
        //unsafe { DoubleField::FCOUNT += 1; }
        self.d -= o.d;
    }

    fn m(&self, o: &DoubleField) -> DoubleField {
        //unsafe { DoubleField::FCOUNT += 1; }
        DoubleField::new(self.d * o.d)
    }

    fn me(&mut self, o: &DoubleField) {
        //unsafe { DoubleField::FCOUNT += 1; }
        self.d *= o.d;
    }

    fn d(&self, o: &DoubleField) -> DoubleField {
        //unsafe { DoubleField::FCOUNT += 1; }
        DoubleField::new(self.d / o.d)
    }

    fn de(&mut self, o: &DoubleField) {
        //unsafe { DoubleField::FCOUNT += 1; }
        self.d /= o.d;
    }

    fn coerce_to_f64(&self) -> f64 {
        self.d
    }
    fn coerce(&self, value: f64) -> DoubleField {
        DoubleField::new(value)
    }

    fn is_zero(&self) -> bool {
        self.d == 0.0
    }

    fn is_one(&self) -> bool {
        self.d == 1.0
    }

    fn zero(&self) -> DoubleField {
        DoubleField::new(0.0)
    }

    fn one(&self) -> DoubleField {
        DoubleField::new(1.0)
    }

    fn copy(&self) -> DoubleField {
        DoubleField::new(self.d)
    }
}


impl IMath for DoubleField {
    fn abs(&self) -> f64 {
        self.d.abs()
    }
    fn sqrt(&mut self) {
        self.d = self.d.sqrt();
    }
}

impl IOrdered for DoubleField {
    fn lt(&self, o: &DoubleField) -> bool {
        self.d < o.d
    }

    fn le(&self, o: &DoubleField) -> bool {
        self.d <= o.d
    }

    fn gt(&self, o: &DoubleField) -> bool {
        self.d > o.d
    }

    fn ge(&self, o: &DoubleField) -> bool {
        self.d >= o.d
    }

    fn e(&self, o: &DoubleField) -> bool {
        self.d == o.d
    }
}

impl ITrigonometric for DoubleField {
    fn sin(&mut self) {
        self.d = self.d.sin();
    }

    fn cos(&mut self) {
        self.d = self.d.cos();
    }
}

impl fmt::Display for DoubleField {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if self.print_short {
            write!(f, "{:6.2}", self.d)
        } else {
            write!(f, "{}", self.d)
        }
    }
}

impl Clone for DoubleField {
    fn clone(&self) -> Self {
        DoubleField::new(self.d)
    }
}

impl Hash for DoubleField {
    fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
        let bits: u64 = self.d.to_bits();
        bits.hash(state);
    }
}

impl PartialEq for DoubleField {
    fn eq(&self, other: &Self) -> bool {
        self.d == other.d
    }
}

impl Eq for DoubleField {}