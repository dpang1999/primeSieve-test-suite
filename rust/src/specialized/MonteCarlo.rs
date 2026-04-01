use crate::helpers::lcg::Lcg;

fn integrate(num_samples: usize) -> f64 {
    let mut rng = Lcg::new(12345, 1345, 16645, 1013904);
    let mut under_curve = 0;
    for _ in 0..num_samples {
        let x: f64 = rng.next_double();
        let y: f64 = rng.next_double();
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
    println!("Rust specialized f64 montecarlo");
    println!("Pi is approximately: {}", pi);
    println!("Num samples: {}", num_samples);
    println!("RMS Error: {}", (std::f64::consts::PI - pi).abs());
}
#[allow(dead_code)]
pub fn run_algorithm() {
    main();
}
