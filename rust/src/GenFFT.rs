use core::fmt;
use rand::Rng;
use crate::generic::complex_field::ComplexField;
use crate::generic::i_field::IField;
use crate::generic::i_primitive_roots::IPrimitiveRoots;
use crate::generic::i_trigonometric::ITrigonometric;
use crate::generic::i_math::IMath;
use crate::generic::i_ordered::IOrdered;
pub mod generic;
use crate::generic::double_field::DoubleField;
use crate::generic::int_mod_p::IntModP;

pub struct GenFFT<N>
where
    N: IField + IMath + IOrdered + IPrimitiveRoots<N> + Clone + fmt::Display,
{
    c: N,
}

impl<N> GenFFT<N>
where
    N: IField + IMath + IOrdered + IPrimitiveRoots<N> + Clone + fmt::Display,
{
    pub fn new(data: N) -> Self {
        Self {
            c: data,
        }
    }

    pub fn transform(&self, data: &mut [N]) {
        self.transform_internal(data, -1);
    }

    pub fn inverse(&self, data: &mut [N]) {
        self.transform_internal(data, 1);
        let nd = data.len();
        let norm = self.c.coerce(nd as f64);

        for d in data.iter_mut() {
            d.de(&norm);
        }
    }

    pub fn test(&self, data: &mut [N]) -> f64 {
        let nd = data.len();
        let mut copy: Vec<N> = data.iter().map(|x| x.copy()).collect();

        self.transform(data);

        println!("After transform: {}", data.iter().map(|x| x.to_string()).collect::<Vec<_>>().join(", "));
        self.inverse(data);
        println!("After inverse: {}", data.iter().map(|x| x.to_string()).collect::<Vec<_>>().join(", "));

      
        let mut diff = 0.0;
        for i in 0..nd {
            let d = data[i].copy();
            let real = d.coerce_to_f64();
            let imag = d.coerce_to_f64();
            let realDiff = real - copy[i].coerce_to_f64();
            let imagDiff = imag - copy[i].coerce_to_f64();
            diff += realDiff * realDiff + imagDiff * imagDiff;
        }
        (diff / (nd*2) as f64).sqrt()

    }

    pub fn make_random(&self, n: usize) -> Vec<ComplexField<N>> {
        let mut data = Vec::with_capacity(n);
        for _ in 0..n {
            data.push(ComplexField::new(
                self.c.coerce(rand::random::<f64>()),
                self.c.coerce(rand::random::<f64>()),
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

    /*fn precompute_roots_of_unity(&self, n: usize, direction: i32) -> Vec<ComplexField<N>> {
        let mut roots = Vec::with_capacity(n);
        for k in 0..n {
            let angle = direction as f64 * 2.0 * std::f64::consts::PI * (k as f64 / n as f64);
            roots.push(ComplexField::new(
                self.c.coerce(angle.cos()),
                self.c.coerce(angle.sin()),
            ));
        }
        roots
    }*/

    fn transform_internal(&self, data: &mut [N], direction: i32) {
        let n = data.len();
        if n == 0 || n == 1 {
            return;
        }
        let logn = Self::log2(n);

        Self::bitreverse(data);
        //println!("After bitreverse: {}", data.iter().map(|x| x.to_string()).collect::<Vec<_>>().join(", "));

        // Precompute roots of unity
        let roots = self.c.precomputeRootsOfUnity(n as u64, direction as u64);

        let mut dual = 1;
        for _bit in 0..logn {
            for a in 0..dual {
                let w = roots[a * (n / (2 * dual))].clone(); // Use precomputed root
                for b in (0..n).step_by(2 * dual) {
                    let i = b + a;
                    let j = b + a + dual;

                    let wd = w.m(&data[j]); // Twiddle factor multiplication
                    data[j] = data[i].s(&wd); // Subtract
                    data[i] = data[i].a(&wd); // Add
                }
            }
            dual *= 2;
            //println!("{}", data[bit]);
        }
    }

    fn bitreverse(data: &mut [N]) {
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
    let fft = GenFFT::new(c);
    let n = 1024;
    
    //let mut data1;
    //let mut data2;
    
    

    let in1 = [38, 0, 44, 87, 6, 45, 22, 93, 0, 0, 0, 0, 0, 0, 0, 0];
    let in2 = [80, 18, 62, 90, 17, 96, 27, 97, 0, 0, 0, 0, 0, 0, 0, 0];
    let out = [3040, 684, 5876, 11172, 5420, 16710, 12546, 20555, 16730, 15704, 21665, 5490, 13887, 4645, 9021, 0];
    let prime = 40961;

    let finite = IntModP::new(0, prime);
    let finiteFFT = GenFFT::new(finite);
    let mut data1 = Vec::with_capacity(in1.len());
    let mut data2 = Vec::with_capacity(in2.len());
    for i in 0..in1.len() {
        data1.push(IntModP::new(in1[i], prime));
        data2.push(IntModP::new(in2[i], prime));
    }
    finiteFFT.transform(&mut data1);
    finiteFFT.transform(&mut data2);

    // Print data1 and data2 after transformation
    println!("data1: {}", data1.iter().map(|x| x.to_string()).collect::<Vec<_>>().join(", "));
    println!("data2: {}", data2.iter().map(|x| x.to_string()).collect::<Vec<_>>().join(", "));

    let mut product = Vec::with_capacity(data1.len());
    for i in 0..data1.len() {
        product.push(data1[i].m(&data2[i]));
    }
    // Print product
    println!("product: {}", product.iter().map(|x| x.to_string()).collect::<Vec<_>>().join(", "));

    // Invert and print product
    finiteFFT.inverse(&mut product);
    println!("product (after inverse): {}", product.iter().map(|x| x.to_string()).collect::<Vec<_>>().join(", "));

    //println!("{}", data.iter().map(|x| x.to_string()).collect::<Vec<_>>().join(", "));
    //println!("n={} => RMS Error={}", n, fft.test(&mut data));
}