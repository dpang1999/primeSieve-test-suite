
use seeded_random::{Random,Seed};

use rust::helpers::lcg::Lcg;
use rust::helpers::prime_sieve::prime_sieve;
use rust::helpers::find_prime::find_prime_congruent_one_mod_n;

pub struct FFT {}


impl FFT
{
    pub fn new() -> Self {
        Self {}
    }

    pub fn transform(&self, data: &mut [f64]) {
        Self::transform_internal(data, -1);
    }

    pub fn inverse(&self, data: &mut [f64]) {
        Self::transform_internal(data, 1);
        let nd = data.len();
        let n = nd/2;
        let norm = 1.0 / n as f64;
        for d in 0..nd {
            data[d] *= norm;
        }
    }

    pub fn test(&self, data: &mut [f64]) -> f64 {
        let nd = data.len();
        let copy: Vec<f64> = data.iter().map(|x| *x).collect();

        self.transform(data);

        println!("After transform: {}", data.iter().map(|x| x.to_string()).collect::<Vec<_>>().join(", "));
        self.inverse(data);
        println!("After inverse: {}", data.iter().map(|x| x.to_string()).collect::<Vec<_>>().join(", "));

        let mut diff = 0.0;
        for i in 0..nd {
            let d = data[i] - copy[i];
            diff += d * d;
        }
        (diff / nd as f64).sqrt()
    }

    pub fn make_random(&self, n: usize) -> Vec<f64> {
        let seed = Seed::unsafe_new(12345);
        let rng = Random::from_seed(seed);
        let nd = 2*n;
        let mut data = Vec::with_capacity(nd);
        for _ in 0..nd {
            data.push(rng.gen());
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

    fn transform_internal(data: &mut [f64], direction: i32) {
        if data.is_empty() {
            return;
        }

        let n = data.len() / 2; // Number of complex elements
        if n == 1 {
            return; // Identity operation
        }

        let logn = Self::log2(n);

        // Bit-reverse the input data for decimation-in-time algorithm
        Self::bitreverse(data);
        //println!("After bit-reverse: {}", data.iter().map(|x| x.to_string()).collect::<Vec<_>>().join(", "));

        // Apply FFT recursion
        for bit in 0..logn {
            let dual = 1 << bit; // 2^bit
            let mut w_real = 1.0;
            let mut w_imag = 0.0;

            let theta = 2.0 * direction as f64 * std::f64::consts::PI / (2.0 * dual as f64);
            let s = theta.sin();
            let t = (theta / 2.0).sin();
            let s2 = 2.0 * t * t;

            // a = 0
            for b in (0..n).step_by(2 * dual) {
                let i = 2 * b;
                let j = 2 * (b + dual);

                let wd_real = data[j];
                let wd_imag = data[j + 1];

                data[j] = data[i] - wd_real;
                data[j + 1] = data[i + 1] - wd_imag;
                data[i] += wd_real;
                data[i + 1] += wd_imag;
            }

            // a = 1 .. (dual-1)
            for a in 1..dual {
                // Trigonometric recurrence for w -> exp(i * theta) * w
                let tmp_real = w_real - s * w_imag - s2 * w_real;
                let tmp_imag = w_imag + s * w_real - s2 * w_imag;
                w_real = tmp_real;
                w_imag = tmp_imag;

                for b in (0..n).step_by(2 * dual) {
                    let i = 2 * (b + a);
                    let j = 2 * (b + a + dual);

                    let z1_real = data[j];
                    let z1_imag = data[j + 1];

                    let wd_real = w_real * z1_real - w_imag * z1_imag;
                    let wd_imag = w_real * z1_imag + w_imag * z1_real;

                    data[j] = data[i] - wd_real;
                    data[j + 1] = data[i + 1] - wd_imag;
                    data[i] += wd_real;
                    data[i + 1] += wd_imag;
                }
            }
            //println!("{} {}", data[bit], data[bit + 1]);
        }
    }

    fn bitreverse(data: &mut [f64]) {
        let n = data.len() / 2; // Number of complex elements
        let nm1 = n - 1;
        let mut i = 0;
        let mut j = 0;
        while i < nm1 {
            if i < j {
                let ii = i * 2; // Real part index
                let jj = j * 2; // Imaginary part index
                data.swap(ii, jj);
                data.swap(ii + 1, jj + 1);
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
    // let mode = 0 be for testing
    let mode = 0;
    let fft = FFT::new();
    if mode != 0 {
        let args: Vec<String> = std::env::args().collect();
        let n: usize = args.get(1).and_then(|s| s.parse().ok()).unwrap_or(16);
        let mut rand = Lcg::new(12345,1345,65,17);
        let mut data: Vec<f64> = Vec::with_capacity(2*n);
        for _ in 0..n {
            let r = rand.next_double();
            data.push(r);
            let i = rand.next_double();
            data.push(i);
        }
        println!("n={} => RMS Error={}", n, fft.test(&mut data));
    }
    else {
       
        let n = 1024;
        let mut data = fft.make_random(n);
        let mut data2: Vec<f64> = vec![
            0.3618031071604718,
            0.932993485288541,
            0.8330913489710237,
            0.32647575623792624,
            0.2355237906476252,
            0.34911535662488336,
            0.4480776326931518,
            0.6381529437838686,
        ];

        //println!("{}", data.iter().map(|x| x.to_string()).collect::<Vec<_>>().join(", "));
        println!("n={} => RMS Error={}", n, fft.test(&mut data));
    }
}