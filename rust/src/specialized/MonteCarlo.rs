use seeded_random::{Random,Seed};

const SEED: u64 = 113;


fn num_flops(num_samples: usize) -> f64 {
    // 3 flops in x^2+y^2 and 1 flop in random routine
    (num_samples as f64) * 4.0
}

fn integrate(num_samples: usize) -> f64 {
    let seed1 = Seed::unsafe_new(SEED);
    let rng = Random::from_seed(seed1);
    let mut under_curve = 0;
    for _ in 0..num_samples {
        let x: f64 = rng.gen();
        let y: f64 = rng.gen();
        if x * x + y * y <= 1.0 {
            under_curve += 1;
        }
    }
    (under_curve as f64 / num_samples as f64) * 4.0
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let mut num_samples = 1_000_000;
    if args.len() > 1 {
        num_samples = args[1].parse().unwrap_or(num_samples);
    }

    let pi = integrate(num_samples);
    println!("Pi is approximately: {}", pi);
    println!("Num samples: {}", num_samples);
    println!("Num flops: {}", num_flops(num_samples));
    println!("RMS Error: {}", (std::f64::consts::PI - pi).abs());
}