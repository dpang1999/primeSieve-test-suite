use std::fmt;
use crate::generic::i_ring::IRing;
use crate::generic::i_invertible::IInvertible;
use crate::generic::i_math::IMath;
use crate::generic::i_ordered::IOrdered;
use crate::generic::i_trigonometric::ITrigonometric;
pub struct DoubleRing {
    pub d: f64,
    pub print_short: bool,
}

impl DoubleRing {
    //pub static mut FCOUNT: i32 = 0;

    pub fn new(d: f64) -> Self {
        DoubleRing { d, print_short: true }
    }

    pub fn copy(&self) -> DoubleRing {
        DoubleRing::new(self.d)
    }

    /* pub fn new_array(size: usize) -> Vec<DoubleRing> {
        vec![DoubleRing::new(0.0); size]
    } */
}


impl IRing for DoubleRing {
    fn a(&self, o: &DoubleRing) -> DoubleRing {
        //unsafe { DoubleRing::FCOUNT += 1; }
        DoubleRing::new(self.d + o.d)
    }

    fn ae(&mut self, o: &DoubleRing) {
        //unsafe { DoubleRing::FCOUNT += 1; }
        self.d += o.d;

    }

    fn s(&self, o: &DoubleRing) -> DoubleRing {
        //unsafe { DoubleRing::FCOUNT += 1; }
        DoubleRing::new(self.d - o.d)
    }

    fn se(&mut self, o: &DoubleRing) {
        //unsafe { DoubleRing::FCOUNT += 1; }
        self.d -= o.d;
    }

    fn m(&self, o: &DoubleRing) -> DoubleRing {
        //unsafe { DoubleRing::FCOUNT += 1; }
        DoubleRing::new(self.d * o.d)
    }

    fn me(&mut self, o: &DoubleRing) {
        //unsafe { DoubleRing::FCOUNT += 1; }
        self.d *= o.d;
    }

    fn coerce_from_i32(&self, i: i32) -> DoubleRing {
        DoubleRing::new(i as f64)
    }

    fn coerce_from_f64(&self, d: f64) -> DoubleRing {
        DoubleRing::new(d)
    }

    fn coerce(&self) -> f64 {
        self.d
    }

    fn is_zero(&self) -> bool {
        self.d == 0.0
    }

    fn is_one(&self) -> bool {
        self.d == 1.0
    }

    fn zero() -> DoubleRing {
        DoubleRing::new(0.0)
    }

    fn one() -> DoubleRing {
        DoubleRing::new(1.0)
    }

    fn copy(&self) -> DoubleRing {
        DoubleRing::new(self.d)
    }
}

impl IInvertible for DoubleRing {
    fn invert(&self) -> DoubleRing {
        //unsafe { DoubleRing::FCOUNT += 1; }
        DoubleRing::new(1.0 / self.d)
    }
}

impl IMath for DoubleRing {
    fn abs(&self) -> DoubleRing {
        DoubleRing::new(self.d.abs())
    }
    fn sqrt(&mut self) {
        self.d = self.d.sqrt();
    }
}

impl IOrdered for DoubleRing {
    fn lt(&self, o: &DoubleRing) -> bool {
        self.d < o.d
    }

    fn le(&self, o: &DoubleRing) -> bool {
        self.d <= o.d
    }

    fn gt(&self, o: &DoubleRing) -> bool {
        self.d > o.d
    }

    fn ge(&self, o: &DoubleRing) -> bool {
        self.d >= o.d
    }
}

impl ITrigonometric for DoubleRing {
    fn sin(&mut self) {
        self.d = self.d.sin();
    }

    fn cos(&mut self) {
        self.d = self.d.cos();
    }
}

impl fmt::Display for DoubleRing {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if self.print_short {
            write!(f, "{:6.2}", self.d)
        } else {
            write!(f, "{}", self.d)
        }
    }
}