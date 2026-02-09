use core::fmt;
use rust::helpers::lcg::Lcg;
use rust::helpers::prime_sieve::prime_sieve;
use rust::helpers::find_prime::find_prime_congruent_one_mod_n;

use crate::generic::complex_field::ComplexField;
use crate::generic::i_field::IField;
use crate::generic::i_primitive_roots::IPrimitiveRoots;
use crate::generic::i_math::IMath;
use crate::generic::i_ordered::IOrdered;
pub mod generic;
use crate::generic::double_field::DoubleField;
use crate::generic::int_mod_p::IntModP;
use crate::generic::int_mod_p::MODULUS;

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
        let copy: Vec<N> = data.iter().map(|x| x.copy()).collect();

        self.transform(data);

        println!("After transform: {}", data.iter().map(|x| x.to_string()).collect::<Vec<_>>().join(", "));
        self.inverse(data);
        println!("After inverse: {}", data.iter().map(|x| x.to_string()).collect::<Vec<_>>().join(", "));

      
        let mut diff = 0.0;
        for i in 0..nd {
            let d = data[i].copy();
            let real = d.coerce_to_f64();
            let imag = d.coerce_to_f64();
            let real_diff = real - copy[i].coerce_to_f64();
            let imag_diff = imag - copy[i].coerce_to_f64();
            diff += real_diff * real_diff + imag_diff * imag_diff;
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
        let roots = self.c.precomputeRootsOfUnity(n as u32, direction);
        //println!("Roots: {}", roots.iter().map(|x| x.to_string()).collect::<Vec<_>>().join(", "));

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
    // let mode = 0 be for testing
    let mode = 1;
    if mode != 0 { 
        // arg 1 = size (N = power of 2)
        // arg 2 = field type (0 = finite field, 1 = complex field)
        let args: Vec<String> = std::env::args().collect();
        let n: usize = args.get(1).and_then(|s| s.parse().ok()).unwrap_or(16);
        let field_type: usize = args.get(2).and_then(|s| s.parse().ok()).unwrap_or(0);
        let mut rand = Lcg::new(12345,1345,65,17);
        if field_type == 0 {
            let mut data1 = Vec::<IntModP>::with_capacity(n);
            //let mut data2 = Vec::<IntModP>::with_capacity(n);
    
            // prime is the largest prime less than num
            let prime = find_prime_congruent_one_mod_n(n);
            MODULUS.set(prime as u128).unwrap();
            for _ in 0..n {
                data1.push(IntModP::new(rand.next_int() as u128 % prime as u128));
                //data2.push(IntModP::new(rand.next_int() as u128 % prime as u128, prime as u128));
            }
            //let data1_clone = data1.clone();
            //println!("Data before transform: {}", data1.iter().map(|x| x.to_string()).collect::<Vec<_>>().join(", "));
            let finite = IntModP::new(0);
            let finite_fft = GenFFT::new(finite);
            println!("Generic Rust FFT Tests");
            println!("Rust Generics, Finite Field, n={}", n);
            for i in 0..10 {
                finite_fft.transform(&mut data1);
                //println!("Data after transform: {}", data1.iter().map(|x| x.to_string()).collect::<Vec<_>>().join(", "));
                finite_fft.inverse(&mut data1);
                //println!("Data after inverse: {}", data1.iter().map(|x| x.to_string()).collect::<Vec<_>>().join(", "));   
                println!("Loop {} done", i);
            }

            /*
            // RMS error
            let mut diff = 0.0;
            for i in 0..n {
                let d = data1[i].copy();
                let val = d.coerce_to_f64();
                let orig = data1_clone[i].coerce_to_f64();
                let real_diff = val - orig;
                diff += real_diff * real_diff;
            }
            let rms = (diff / n as f64).sqrt();
            println!("RMS Error: {}", rms);    */    

        }
        else if field_type == 1 {
            let c = ComplexField::new(DoubleField::new(0.0), DoubleField::new(0.0));
            let fft = GenFFT::new(c);
            let mut data1 = Vec::<ComplexField<DoubleField>>::with_capacity(n);
            for _ in 0..n {
                data1.push(ComplexField::new(
                    DoubleField::new(rand.next_double()),
                    DoubleField::new(rand.next_double()),
                ));
            }

            println!("Generic Rust FFT Tests");
            println!("Rust Generics, Complex Field, n={}", n);
            for i in 0..10 {
                fft.transform_internal(&mut data1, -1);
                fft.transform_internal(&mut data1, 1);
                println!("Loop {} done", i);
            }

            //let data1_clone = data1.clone();
            //println!("Data before transform: {}", data1.iter().map(|x| x.to_string()).collect::<Vec<_>>().join(", "));
            //fft.transform(&mut data1);
            //println!("Data after transform: {}", data1.iter().map(|x| x.to_string()).collect::<Vec<_>>().join(", "));
            //fft.inverse(&mut data1);
            //println!("Data after inverse: {}", data1.iter().map(|x| x.to_string()).collect::<Vec<_>>().join(", "));

            // RMS error
            /*let mut diff = 0.0;
            for i in 0..n {
                let d = data1[i].copy();
                let real = d.re.coerce_to_f64();
                let imag = d.im.coerce_to_f64();
                let real_diff = real - data1_clone[i].re.coerce_to_f64();
                let imag_diff = imag - data1_clone[i].im.coerce_to_f64();
                diff += real_diff * real_diff + imag_diff * imag_diff;
            }
            let rms = (diff / (n*2) as f64).sqrt();
            println!("RMS Error: {}", rms);*/


        }

    }
    else {
        // test: 0 for finite field, 1 complex field
        let test = 3;    
        if test == 0 {
            
            let in1 = [38, 0, 44, 87, 6, 45, 22, 93, 0, 0, 0, 0, 0, 0, 0, 0];
            let in2 = [80, 18, 62, 90, 17, 96, 27, 97, 0, 0, 0, 0, 0, 0, 0, 0];
            //let out = [3040, 684, 5876, 11172, 5420, 16710, 12546, 20555, 16730, 15704, 21665, 5490, 13887, 4645, 9021, 0];
            let prime = 40961;
            MODULUS.set(prime as u128).unwrap();
            
            
            /* 
            let in1: [u128; 16] = [
                11400, 28374, 23152, 9576, 29511, 20787, 13067, 14015, 0, 0, 0, 0, 0, 0, 0, 0,
            ];
            let in2: [u128; 16] = [
                30268, 20788, 8033, 15446, 26275, 11619, 2494, 7016, 0, 0, 0, 0, 0, 0, 0, 0,
            ];
            /*let out: [u128; 16] = [
                345055200, 1095807432, 1382179648, 1175142886, 2016084656, 2555168834,
                2179032777, 1990011337, 1860865174, 1389799087, 942120918, 778961552,
                341270975, 126631482, 98329240, 0
            ];*/
            let prime: u128 = 3221225473;
            */
        /*
            let in1: [u128; 64] = [33243586, 638827078, 767661659, 778933286, 790244973, 910208076, 425757125,
                478004096, 153380495, 205851834, 668901196, 15731080, 899763115, 551605421,
                181279081, 600279047, 711828654, 483031418, 737709105, 20544909, 609397212,
                201989947, 215952988, 206613081, 471852626, 889775274, 992608567, 947438771,
                969970961, 676943009, 934992634, 922939225, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];
            let in2: [u128; 64] = [194132110, 219972873, 66644114, 902841100, 565039275, 540721923, 810650854,
                702680360, 147944788, 859947137, 59055854, 288190067, 537655879, 836782561,
                308822170, 315498953, 417177801, 640439652, 198304612, 525827778, 115633328,
                285831984, 136721026, 203065689, 884961191, 222965182, 735241234, 746745227,
                667772468, 739110962, 610860398, 965331182, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];
            /*let out: [u128; 64] = [6453647494146460, 131329535698517158, 291767894660778388, 392668443347293259,
                971459521481104784, 1474458811520325621, 1844928110064910283, 2357021332184901128,
                2928892267161886295, 2725517850003984528, 3202505799926570519, 2918543444592941968,
                2772488376791744089, 3248633108357294538, 3254615389814072180, 3638020871734883400,
                55160505208503622, 3969469665294621400, 439789777768675993, 916737048670338429,
                157193402339279849, 1030499289809835368, 534708807109284987, 462608833776141716,
                518270737313306417, 990302136704222252, 862673986833243374, 1706781055673683080,
                2148213235654123180, 4027029548560043607, 3715706394243238489, 966330325631268533,
                724857759400778139, 1014165568394318451, 978244158856038395, 3518954508900415555,
                3481727912868647859, 2905676401026905092, 1913454655595000205, 2281030150295966751,
                2048468707271352286, 1955651308030723278, 1936345891479581000, 2116568874488615349,
                1964776204460631657, 594938508019154838, 665031798826217600, 435270820221219547,
                3944115800695200119, 3877068415832542765, 3375534600145876311, 3739051895812367546,
                3787681810231019302, 3846806706428246918, 215267241912496193, 433277273552403593,
                32647322247915044, 4082693161306839314, 3321007834415954245, 2657237599459774692,
                1906778666014199420, 1466364566853824938, 890942012983413950, 0];
            */
            let prime: u128 = 4179340454199820289;*/

            let finite = IntModP::new(0);
            let finite_fft = GenFFT::new(finite);
            let mut data1 = Vec::with_capacity(in1.len());
            let mut data2 = Vec::with_capacity(in2.len());
            for i in 0..in1.len() {
                data1.push(IntModP::new(in1[i]));
                data2.push(IntModP::new(in2[i]));
            }
            finite_fft.transform(&mut data1);
            finite_fft.transform(&mut data2);

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
            finite_fft.inverse(&mut product);
            println!("product (after inverse): {}", product.iter().map(|x| x.to_string()).collect::<Vec<_>>().join(", "));
        }
        else if test == 1 {
            let c = ComplexField::new(DoubleField::new(0.0), DoubleField::new(0.0));
            let fft = GenFFT::new(c);
            let n = 4;

            let mut data1 = Vec::<ComplexField<DoubleField>>::with_capacity(n);
            let mut data2 = Vec::<ComplexField<DoubleField>>::with_capacity(n);

            data1.push(ComplexField::new(DoubleField::new(0.3618031071604718), DoubleField::new(0.932993485288541)));
            data1.push(ComplexField::new(DoubleField::new(0.8330913489710237), DoubleField::new(0.32647575623792624)));
            data1.push(ComplexField::new(DoubleField::new(0.2355237906476252), DoubleField::new(0.34911535662488336)));
            data1.push(ComplexField::new(DoubleField::new(0.4480776326931518), DoubleField::new(0.6381529437838686)));


            println!("{}", data1.iter().map(|x| x.to_string()).collect::<Vec<_>>().join(", "));
            println!("n={} => RMS Error={}", n, fft.test(&mut data1));

            //let in1 = [38, 0, 44, 87, 6, 45, 22, 93, 0, 0, 0, 0, 0, 0, 0, 0];
            //let in2 = [80, 18, 62, 90, 17, 96, 27, 97, 0, 0, 0, 0, 0, 0, 0, 0];
            //let in1 = [11400, 28374, 23152, 9576, 29511, 20787, 13067, 14015, 0, 0, 0, 0, 0, 0, 0, 0];
            //let in2 = [30268, 20788, 8033, 15446, 26275, 11619, 2494, 7016, 0, 0, 0, 0, 0, 0, 0, 0];
            let in1 = [33243586, 638827078, 767661659, 778933286, 790244973, 910208076, 425757125,
                478004096, 153380495, 205851834, 668901196, 15731080, 899763115, 551605421,
                181279081, 600279047, 711828654, 483031418, 737709105, 20544909, 609397212,
                201989947, 215952988, 206613081, 471852626, 889775274, 992608567, 947438771,
                969970961, 676943009, 934992634, 922939225, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];
            let in2 = [194132110, 219972873, 66644114, 902841100, 565039275, 540721923, 810650854,
                702680360, 147944788, 859947137, 59055854, 288190067, 537655879, 836782561,
                308822170, 315498953, 417177801, 640439652, 198304612, 525827778, 115633328,
                285831984, 136721026, 203065689, 884961191, 222965182, 735241234, 746745227,
                667772468, 739110962, 610860398, 965331182, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];


            let mut data1 = Vec::with_capacity(in1.len() / 2);
            let mut data2 = Vec::with_capacity(in2.len() / 2);
            for i in (0..in1.len()).step_by(2) {
                data1.push(ComplexField::new(DoubleField::new(in1[i] as f64), DoubleField::new(in1[i+1] as f64)));
                data2.push(ComplexField::new(DoubleField::new(in2[i] as f64), DoubleField::new(in2[i+1] as f64)));
            }
            // Print data1 and data2 before transformation
            println!("data1 (before transform): {}", data1.iter().map(|x| x.to_string()).collect::<Vec<_>>().join(", "));
            println!("data2 (before transform): {}", data2.iter().map(|x| x.to_string()).collect::<Vec<_>>().join(", "));
            fft.transform(&mut data1);
            fft.transform(&mut data2);

            // Print data1 and data2 after transformation
            println!("data1 (after transform): {}", data1.iter().map(|x| x.to_string()).collect::<Vec<_>>().join(", "));
            println!("data2 (after transform): {}", data2.iter().map(|x| x.to_string()).collect::<Vec<_>>().join(", "));

            let mut product = Vec::with_capacity(data1.len());
            for i in 0..data1.len() {
                product.push(data1[i].m(&data2[i]));
            }
            // Print product
            println!("product: {}", product.iter().map(|x| x.to_string()).collect::<Vec<_>>().join(", "));

            // Invert and print product
            fft.inverse(&mut product);
            println!("product (after inverse): {}", product.iter().map(|x| x.to_string()).collect::<Vec<_>>().join(", "));
        }
        else {
            let mut rand = Lcg::new(12345,1345,65,17);
            let mut random_numbers = [0; 10];
            let mut random_doubles = [0.0; 10];
            for i in 0..10 {
                random_numbers[i] = rand.next_int();
                random_doubles[i] = rand.next_double();  
            }
            println!("Random Integers: {}", random_numbers.iter().map(|x| x.to_string()).collect::<Vec<_>>().join(", "));
            println!("Random Doubles: {}", random_doubles.iter().map(|x| x.to_string()).collect::<Vec<_>>().join(", "));
        }
    }
}