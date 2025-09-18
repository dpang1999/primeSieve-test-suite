use core::fmt;
use std::f64::consts::PI;
use crate::generic::i_field::IField;
use crate::generic::i_ordered::IOrdered;
use crate::generic::i_math::IMath;
use crate::generic::single_field::SingleField;
use crate::generic::double_field::DoubleField;
use crate::generic::int_mod_p::IntModP;
use crate::generic::i_primitive_roots::IPrimitiveRoots;

pub struct ComplexField<T> {
    pub re: T,
    pub im: T,
}

impl<T: IField + IOrdered> ComplexField<T> {
    pub fn new(re: T, im: T) -> Self {
        ComplexField { re, im }
    }
}

impl <T: IField + IOrdered> IField for ComplexField<T> {
    fn copy(&self) -> Self {
        ComplexField::new(self.re.copy(), self.im.copy())
    }

    fn a(&self, o: &ComplexField<T>) -> ComplexField<T> {
        ComplexField::new(self.re.a(&o.re), self.im.a(&o.im))
    }

    fn ae(&mut self, o: &ComplexField<T>) {
        self.re = self.re.a(&o.re);
        self.im = self.im.a(&o.im);
    }

    fn s(&self, o: &ComplexField<T>) -> ComplexField<T> {
        ComplexField::new(self.re.s(&o.re), self.im.s(&o.im))
    }

    fn se(&mut self, o: &ComplexField<T>) {
        self.re = self.re.s(&o.re);
        self.im = self.im.s(&o.im);
    }

    fn m(&self, o: &ComplexField<T>) -> ComplexField<T> {
        // (a + bi) * (c + di) = (ac - bd) + (ad + bc)i
        ComplexField::new(
            self.re.m(&o.re).s(&self.im.m(&o.im)),
            self.re.m(&o.im).a(&self.im.m(&o.re)),
        )
    }

    fn me(&mut self, o: &ComplexField<T>) {
        let temp_re = self.re.m(&o.re).s(&self.im.m(&o.im));
        let temp_im = self.re.m(&o.im).a(&self.im.m(&o.re));
        self.re = temp_re;
        self.im = temp_im;
    }

    fn d(&self, o: &ComplexField<T>) -> ComplexField<T> {
        // (a + bi) / (c + di) = [(ac + bd) / (c^2 + d^2)] + [(bc - ad) / (c^2 + d^2)]i
        let denom = o.re.m(&o.re).a(&o.im.m(&o.im));
        ComplexField::new(
            self.re.m(&o.re).a(&self.im.m(&o.im)).d(&denom),
            self.im.m(&o.re).s(&self.re.m(&o.im)).d(&denom),
        )
    }

    fn de(&mut self, o: &ComplexField<T>) {
        let denom = o.re.m(&o.re).a(&o.im.m(&o.im));
        let temp_re = self.re.m(&o.re).a(&self.im.m(&o.im)).d(&denom);
        let temp_im = self.im.m(&o.re).s(&self.re.m(&o.im)).d(&denom);
        self.re = temp_re;
        self.im = temp_im;
    }

    fn coerce_to_f64(&self) -> f64 {
        (self.re.coerce_to_f64().powi(2) + self.im.coerce_to_f64().powi(2)).sqrt()
    }
    fn coerce(&self, value: f64) -> ComplexField<T> {
        ComplexField::new(self.re.coerce(value), self.im.coerce(0.0))
    }

    fn is_zero(&self) -> bool {
        self.re.is_zero() && self.im.is_zero()
    }
    fn is_one(&self) -> bool {
        self.re.is_one() && self.im.is_zero()
    }
    fn zero(&self) -> ComplexField<T> {
        ComplexField::new(T::zero(&self.re), T::zero(&self.im))
    }
    fn one(&self) -> ComplexField<T> {
        ComplexField::new(T::one(&self.re), T::zero(&self.im))
    }
}

impl<T: IField + IOrdered> IOrdered for ComplexField<T> {
    fn lt(&self, o: &ComplexField<T>) -> bool {
        self.re.lt(&o.re) || (self.re.eq(&o.re) && self.im.lt(&o.im))
    }

    fn le(&self, o: &ComplexField<T>) -> bool {
        self.re.lt(&o.re) || (self.re.eq(&o.re) && (self.im.lt(&o.im) || self.im.eq(&o.im)))
    }

    fn gt(&self, o: &ComplexField<T>) -> bool {
        self.re.gt(&o.re) || (self.re.eq(&o.re) && self.im.gt(&o.im))
    }

    fn ge(&self, o: &ComplexField<T>) -> bool {
        self.re.gt(&o.re) || (self.re.eq(&o.re) && (self.im.gt(&o.im) || self.im.eq(&o.im)))
    }

    fn eq(&self, o: &ComplexField<T>) -> bool {
        self.re.eq(&o.re) && self.im.eq(&o.im)
    }
}

impl<T: IField + IOrdered + fmt::Display> fmt::Display for ComplexField<T> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if self.im.is_zero() {
            write!(f, "{}", self.re)
        } else if self.re.is_zero() {
            write!(f, "{}i", self.im)
        } else {
            write!(
                f,
                "{}{}{}i",
                self.re,
                if self.im.lt(&T::zero(&self.im)) { "" } else { "+" },
                self.im
            )
        }
    }
}

impl<T: IField + IOrdered + IMath> IMath for ComplexField<T> {
    fn abs(&self) -> f64 {
        let re = self.re.copy();
        let im = self.im.copy();
        let mut temp = re.m(&re).a(&im.m(&im));
        temp.sqrt();
        temp.abs()
    }

    fn sqrt(&mut self) {
        // Square root of a complex number is not straightforward and is not implemented here.
        panic!("Square root not implemented for ComplexField");
    }
}

impl<T> Clone for ComplexField<T>
where
    T: IField + IOrdered + Clone,
{
    fn clone(&self) -> Self {
        ComplexField::new(self.re.clone(), self.im.clone())
    }
}



// Implement PrimitiveRoot for SingleField numbers
impl IPrimitiveRoots<ComplexField<SingleField>> for ComplexField<SingleField> {
    fn primitive_root(&self, n: u128) -> Self {
        if n == 0 {
            panic!("n must be positive");
        }

        // Compute the angle for the primitive root of unity
        let angle = 2.0 * PI / n as f64;

        // Compute real and imaginary parts
        let real = SingleField::new(angle.cos() as f32);
        let imag = SingleField::new(angle.sin() as f32);

        ComplexField::new(real, imag)
    }

    fn pow(&self, exponent: u128) -> Self {
        if exponent == 0 {
            return self.one(); // Any number to the power of 0 is 1
        }

        // Convert to polar form
        let r = self.abs(); // Modulus
        let theta = self.im.coerce_to_f64().atan2(self.re.coerce_to_f64()); // Argument (angle)

        // Compute new modulus and argument
        let new_r = r.powi(exponent as i32); // r^exponent
        let new_theta = theta * exponent as f64; // theta * exponent

        // Convert back to rectangular form
        let real = self.re.coerce(new_r * new_theta.cos());
        let imag = self.im.coerce(new_r * new_theta.sin());
        ComplexField::new(real, imag)
    }

    fn precomputeRootsOfUnity(&self, n: u32, direction: i32) -> Vec<ComplexField<SingleField>> {
        let mut roots = Vec::new();
        for k in 0..n {
            let angle = 2.0 * PI * k as f64 / n as f64 * direction as f64;
            let real = SingleField::new(angle.cos() as f32);
            let imag = SingleField::new(angle.sin() as f32);
            roots.push(ComplexField::new(real, imag));
        }
        roots
    }
}

// Implement PrimitiveRoot for DoubleField numbers
impl IPrimitiveRoots<ComplexField<DoubleField>> for ComplexField<DoubleField> {
    fn primitive_root(&self, n: u128) -> Self {
        if n == 0 {
            panic!("n must be positive");
        }

        // Compute the angle for the primitive root of unity
        let angle = 2.0 * PI / n as f64;

        // Compute real and imaginary parts
        let real = DoubleField::new(angle.cos());
        let imag = DoubleField::new(angle.sin());

        ComplexField::new(real, imag)
    }

    //Polar form exponentiation
    fn pow(&self, exponent: u128) -> Self {
        if exponent == 0 {
            return self.one(); // Any number to the power of 0 is 1
        }


        // Convert to polar form
        let r = self.abs(); // Modulus
        let theta = self.im.coerce_to_f64().atan2(self.re.coerce_to_f64()); // Argument (angle)

        // Compute new modulus and argument
        let new_r = r.powi(exponent as i32); // r^exponent
        let new_theta = theta * exponent as f64; // theta * exponent

        // Convert back to rectangular form
        let real = self.re.coerce(new_r * new_theta.cos());
        let imag = self.im.coerce(new_r * new_theta.sin());
        ComplexField::new(real, imag)
    }

    fn precomputeRootsOfUnity(&self, n: u32, direction: i32) -> Vec<ComplexField<DoubleField>> {
        let mut roots = Vec::new();
        for k in 0..n {
            let angle = 2.0 * PI * k as f64 / n as f64 * direction as f64;
            let real = self.re.coerce(angle.cos());
            let imag = self.im.coerce(angle.sin());
            roots.push(ComplexField::new(real, imag));
        }
        roots
    }
}

// Implement PrimitiveRoot for finite fields
impl IPrimitiveRoots<ComplexField<IntModP>> for ComplexField<IntModP> {
    fn primitive_root(&self, n: u128) -> Self {
        if n == 0 {
            panic!("n must be positive");
        }

        // For finite fields, compute the primitive root algebraically
        let real = self.re.primitive_root(n); // Assuming T has a primitive_root method
        let imag = self.im.zero();

        ComplexField::new(real, imag)
    }

    fn pow(&self, exponent: u128) -> Self {
        if exponent == 0 {
            return self.one(); // Any number to the power of 0 is 1
        }

        // Convert to polar form
        let r = self.abs(); // Modulus
        let theta = self.im.coerce_to_f64().atan2(self.re.coerce_to_f64()); // Argument (angle)

        // Compute new modulus and argument
        let new_r = r.powi(exponent as i32); // r^exponent
        let new_theta = theta * exponent as f64; // theta * exponent

        // Convert back to rectangular form
        let real = self.re.coerce(new_r * new_theta.cos());
        let imag = self.im.coerce(new_r * new_theta.sin());
        ComplexField::new(real, imag)
    }

    fn precomputeRootsOfUnity(&self, n: u32, direction: i32) -> Vec<ComplexField<IntModP>> {
        let mut roots = Vec::new();
        for k in 0..n {
            let angle = 2.0 * PI * k as f64 / n as f64 * direction as f64;
            let real = self.re.coerce(angle.cos());
            let imag = self.im.coerce(angle.sin());
            roots.push(ComplexField::new(real, imag));
        }
        roots
    }
}