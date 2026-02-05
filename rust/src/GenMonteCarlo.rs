use seeded_random::{Random,Seed};
use crate::generic::double_field::DoubleField;
use crate::generic::single_field::SingleField;
use crate::generic::int_mod_p::IntModP;
use crate::generic::int_mod_p::MODULUS;
use crate::generic::i_field::IField;
use crate::generic::i_ordered::IOrdered;
pub mod generic;
const SEED: u64 = 113;


fn num_flops(num_samples: usize) -> f64 {
    // 3 flops in x^2+y^2 and 1 flop in random routine
    (num_samples as f64) * 4.0
}

fn integrate<T: IField + IOrdered>(t: &T, num_samples: usize) -> f64 {
    let seed1 = Seed::unsafe_new(SEED);
    let rng = Random::from_seed(seed1);
    let mut under_curve = 0;
    for _ in 0..num_samples {
        let x = t.coerce(rng.gen());
        let y = t.coerce(rng.gen());
        if x.m(&x).a(&y.m(&y)).le(&t.one()) {
            under_curve += 1;
        }
    }
    (under_curve as f64 / num_samples as f64) * 4.0
}

fn main() {
    // arg1 = num_samples
    // arg2 = mode (1=SingleField, 2=DoubleField, else IntModP)
    let args: Vec<String> = std::env::args().collect();
    let mut num_samples = 1_000_000;
    let mut mode = 0; // 1 for SingleField, else for DoubleField. Don't think IntModP makes much sense here
    let prime = 1_000_000_007;
    MODULUS.set(prime).unwrap();
    if args.len() > 1 {
        num_samples = args[1].parse().unwrap_or(1000000);
    }
    if args.len() > 2 {
        mode = args[2].parse().unwrap_or(2);
    }
    let pi;
    if(mode == 1) {
        let temp = SingleField::new(0.0);
        pi = integrate(&temp, num_samples);
    }
    else if (mode == 2) {
        let temp = DoubleField::new(0.0);
        pi = integrate(&temp, num_samples);
    }
    else {
        let temp = IntModP::new(0);
        pi = integrate(&temp, num_samples);
    }
    println!("Pi is approximately: {}", pi);
    println!("Num samples: {}", num_samples);
    println!("Num flops: {}", num_flops(num_samples));
    println!("RMS Error: {}", (std::f64::consts::PI - pi).abs());
}