use core::fmt;
use rand::Rng;
use crate::generic::complex_field::ComplexField;
use crate::generic::i_field::IField;
use crate::generic::i_trigonometric::ITrigonometric;
use crate::generic::i_math::IMath;
use crate::generic::i_ordered::IOrdered;
pub mod generic;
use crate::generic::double_field::DoubleField;

pub struct GenFFT<N>
where
    N: IField + ITrigonometric + IMath + IOrdered + Clone + fmt::Display,
{
    c: ComplexField<N>,
}

impl<N> GenFFT<N>
where
    N: IField + ITrigonometric + IMath + IOrdered + Clone + fmt::Display,
{
    pub fn new(re: N, im: N) -> Self {
        Self {
            c: ComplexField::new(re, im),
        }
    }

    pub fn transform(&self, data: &mut [ComplexField<N>]) {
        self.transform_internal(data, -1);
    }

    pub fn inverse(&self, data: &mut [ComplexField<N>]) {
        self.transform_internal(data, 1);
        let nd = data.len();
        let norm = self.c.coerce(1.0 / nd as f64);

        for d in data.iter_mut() {
            d.me(&norm);
        }
    }

    pub fn test(&self, data: &mut [ComplexField<N>]) -> f64 {
        let nd = data.len();
        let mut copy: Vec<ComplexField<N>> = data.iter().map(|x| x.copy()).collect();

        self.transform(data);

        //println!("After transform: {}", data.iter().map(|x| x.to_string()).collect::<Vec<_>>().join(", "));
        self.inverse(data);
        //println!("After inverse: {}", data.iter().map(|x| x.to_string()).collect::<Vec<_>>().join(", "));

      
        let mut diff = 0.0;
        for i in 0..nd {
            let d = data[i].copy();
            let real = d.re.coerce_to_f64();
            let imag = d.im.coerce_to_f64();
            let realDiff = real - copy[i].re.coerce_to_f64();
            let imagDiff = imag - copy[i].im.coerce_to_f64();
            diff += realDiff * realDiff + imagDiff * imagDiff;
        }
        (diff / (nd*2) as f64).sqrt()

    }

    pub fn make_random(&self, n: usize) -> Vec<ComplexField<N>> {
        let mut data = Vec::with_capacity(n);
        for _ in 0..n {
            data.push(ComplexField::new(
                self.c.re.coerce(rand::random::<f64>()),
                self.c.re.coerce(rand::random::<f64>()),
            ));
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

    fn transform_internal(&self, data: &mut [ComplexField<N>], direction: i32) {
        let n = data.len();
        if n == 0 || n == 1 {
            return;
        }
        let logn = Self::log2(n);

        Self::bitreverse(data);
        //println!("After bitreverse: {}", data.iter().map(|x| x.to_string()).collect::<Vec<_>>().join(", "));

        let n1 = self.c.re.coerce(1.0);
        let n2 = self.c.re.coerce(2.0);
        let c10 = ComplexField::new(self.c.re.coerce(1.0), self.c.re.coerce(0.0));
        let mut dual = 1;
        for bit in 0..logn {
            let mut w = c10.copy();
            let mut theta = self.c.re.coerce(
                2.0 * direction as f64 * std::f64::consts::PI / (2.0 * dual as f64),
            );
            let mut s = theta.copy();
            s.sin();
            let mut t = theta.copy();
            t.de(&n2);
            t.sin();
            
            let mut s2 = t.copy();
            let temp = s2.copy();
            s2.me(&temp);
            s2.me(&n2);

            for b in (0..n).step_by(2 * dual) {
                let i = b;
                let j = b + dual;
                let wd = data[j].copy();
                let mut tmp = data[i].copy();
                tmp.se(&wd);
                data[j] = tmp;
                data[i].ae(&wd);
            }

            for a in 1..dual {
                w.me(&ComplexField::new(n1.s(&s2), s.clone()));
                for b in (0..n).step_by(2 * dual) {
                    let i = b + a;
                    let j = b + a + dual;
                    let wd = w.m(&data[j]);
                    data[j] = data[i].s(&wd);
                    data[i].ae(&wd);
                }
            }
            dual *= 2;
            println!("{}", data[bit]);
        }
    }

    fn bitreverse(data: &mut [ComplexField<N>]) {
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
    let c = ComplexField::new(DoubleField::new(0.0), DoubleField::new(0.0));
    let fft = GenFFT::new(DoubleField::new(0.0), DoubleField::new(0.0));
    let n = 4;
    let mut data = fft.make_random(n);
    let mut data2 = Vec::with_capacity(4);
    for i in 0..4 {
        data2.push(ComplexField::new(
            DoubleField::new(2 as f64),
            DoubleField::new(0 as f64),
        ));
    }
    let mut data3 = Vec::with_capacity(4);
    data3.push(ComplexField::new(
        DoubleField::new(0.3618031071604718),
        DoubleField::new(0.932993485288541),
    ));
    data3.push(ComplexField::new(
        DoubleField::new(0.8330913489710237),
        DoubleField::new(0.32647575623792624),
    ));
    data3.push(ComplexField::new(
        DoubleField::new(0.2355237906476252),
        DoubleField::new(0.34911535662488336),
    ));
    data3.push(ComplexField::new(
        DoubleField::new(0.4480776326931518),
        DoubleField::new(0.6381529437838686),
    ));
    println!("{}", data3.iter().map(|x| x.to_string()).collect::<Vec<_>>().join(", "));
    println!("n={} => RMS Error={}", n, fft.test(&mut data3));
}