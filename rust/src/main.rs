pub mod helpers;
fn main() {
    let args: Vec<String> = std::env::args().collect();
    let mut max: usize = 42;
    if args.len() > 1 { // program name is always the first argument
        max = args[1].parse().unwrap();
    }
    
    let temp: Vec<bool> = prime_sieve(max);
    for i in 2..temp.len() {
        if temp[i] {
            println!("{}",i);

        }
    }
    
}

fn prime_sieve(num:usize) -> Vec<bool> {
    let mut numbers:Vec<bool> = vec![true;num];
    numbers[0] = false;
    numbers[1] = false;
    for i in 2..num {
        //interestingly i never takes the value of num, non-inclusive end range
        if numbers[i] {
            let mut j: usize = i;
            let mut current: usize = j*i;
            while current<num {
                numbers[current] = false;
                j+=1;
                current = j*i;
            }
        
        }
    }
    numbers
}