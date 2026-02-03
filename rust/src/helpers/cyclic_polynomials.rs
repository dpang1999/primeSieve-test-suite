fn cyclic_polynomials(n: usize) -> Vec<String> {
    let vars: Vec<String> = (0..=n-1).map(|i| format!("x{}", i)).collect();
    let mut polys = Vec::new();

    for k in 1..n {
        let mut terms = Vec::new();
        for start in 0..n {
            let mut term = String::new();
            for offset in 0..k {
                let idx = (start + offset) % n;
                if offset == 0 {
                    term.push_str(&vars[idx]);
                } else {
                    term.push('*');
                    term.push_str(&vars[idx]);
                }
            }
            terms.push(term);
        }
        polys.push(terms.join(" + "));
    }
    let prod = vars.join("*");
    polys.push(format!("{} - 1", prod));

    polys
}

fn main() {
    let system = cyclic_polynomials(3);
    for (i, poly) in system.iter().enumerate() {
        println!("f{} = {}", i+1, poly);
    }
}