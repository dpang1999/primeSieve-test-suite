use std::fmt;
use crate::generic::i_field::IField;
use crate::generic::i_math::IMath;
use crate::generic::i_ordered::IOrdered;
use crate::generic::i_trigonometric::ITrigonometric;
pub struct SingleField {
    pub f: f32,
    pub print_short: bool,
}

impl SingleField {
    //pub static mut FCOUNT: i32 = 0;

    pub fn new(f: f32) -> Self {
        SingleField { f, print_short: true }
    }

    pub fn copy(&self) -> SingleField {
        SingleField::new(self.f)
    }

    /* pub fn new_array(size: usize) -> Vec<SingleField> {
        vec![SingleField::new(0.0); size]
    } */
}


impl IField for SingleField {
    fn a(&self, o: &SingleField) -> SingleField {
        //unsafe { SingleField::FCOUNT += 1; }
        SingleField::new(self.f + o.f)
    }

    fn ae(&mut self, o: &SingleField) {
        //unsafe { SingleField::FCOUNT += 1; }
        self.f += o.f;

    }

    fn s(&self, o: &SingleField) -> SingleField {
        //unsafe { SingleField::FCOUNT += 1; }
        SingleField::new(self.f - o.f)
    }

    fn se(&mut self, o: &SingleField) {
        //unsafe { SingleField::FCOUNT += 1; }
        self.f -= o.f;
    }

    fn m(&self, o: &SingleField) -> SingleField {
        //unsafe { SingleField::FCOUNT += 1; }
        SingleField::new(self.f * o.f)
    }

    fn me(&mut self, o: &SingleField) {
        //unsafe { SingleField::FCOUNT += 1; }
        self.f *= o.f;
    }

    fn d(&self, o: &SingleField) -> SingleField {
        //unsafe { SingleField::FCOUNT += 1; }
        SingleField::new(self.f / o.f)
    }

    fn de(&mut self, o: &SingleField) {
        //unsafe { SingleField::FCOUNT += 1; }
        self.f /= o.f;
    }

    /*fn coerce(&self) -> f64 {
        self.f as f64
    }*/

    fn is_zero(&self) -> bool {
        self.f == 0.0
    }

    fn is_one(&self) -> bool {
        self.f == 1.0
    }

    fn zero(&self) -> SingleField {
        SingleField::new(0.0)
    }

    fn one(&self) -> SingleField {
        SingleField::new(1.0)
    }

    fn copy(&self) -> SingleField {
        SingleField::new(self.f)
    }
}


impl IMath for SingleField {
    fn abs(&self) -> SingleField {
        SingleField::new(self.f.abs())
    }
    fn sqrt(&mut self) {
        self.f = self.f.sqrt();
    }
}

impl IOrdered for SingleField {
    fn lt(&self, o: &SingleField) -> bool {
        self.f < o.f
    }

    fn le(&self, o: &SingleField) -> bool {
        self.f <= o.f
    }

    fn gt(&self, o: &SingleField) -> bool {
        self.f > o.f
    }

    fn ge(&self, o: &SingleField) -> bool {
        self.f >= o.f
    }

    fn eq(&self, o: &SingleField) -> bool {
        self.f == o.f
    }
}

impl ITrigonometric for SingleField {
    fn sin(&mut self) {
        self.f = self.f.sin();
    }

    fn cos(&mut self) {
        self.f = self.f.cos();
    }
}

impl fmt::Display for SingleField {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if self.print_short {
            write!(f, "{:6.2}", self.f)
        } else {
            write!(f, "{}", self.f)
        }
    }
}