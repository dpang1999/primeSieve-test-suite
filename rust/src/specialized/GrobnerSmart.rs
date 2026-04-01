use std::collections::HashSet;

use std::hash::{Hash, Hasher};

static mut TERM_ORDER: TermOrder = TermOrder::Lex; // default to lex order, can be set to GrLex or RevLex as well
static mut MODULUS: u64 = 7; // default modulus for coefficients, can be changed as needed

use crate::helpers::lcg::Lcg;


#[derive(Clone, Debug, PartialEq)]
pub struct Term {
    pub coefficient: u64,
    pub exponents: u64, // Bitpacking of exponents, 8 bits per variable, 6 variables, last 16 bits for degree
}

#[derive(Copy, Clone, Debug, PartialEq)]
pub enum TermOrder {
    Lex,
    GrLex,
    RevLex
}

impl Term {
    /// Create a new Term with human-readable exponents (array of 6 variables) and degree computed automatically.
    pub fn from_exponents(coefficient: u64, exponents: [u8; 6]) -> Self {
        let mut bitpacked_exponents: u64 = 0;

        // Pack the exponents with the first exponent at the left (bits 47..40)
        // Layout: [63..48]=degree, [47..40]=e0, [39..32]=e1, ... [7..0]=e5
        for (i, &exp) in exponents.iter().enumerate() {
            let shift = 40_u32.saturating_sub((8 * i) as u32);
            bitpacked_exponents |= (exp as u64) << shift;
        }

        // Compute the total degree (sum of all exponents)
        let degree: u16 = exponents.iter().map(|&e| e as u16).sum();

        // Add the degree to the top 16 bits of the u64
        bitpacked_exponents |= (degree as u64) << 48;

        Term {
            coefficient,
            exponents: bitpacked_exponents,
        }
    }
    /* // Helper: pack an array of 6 exponents into the 48-bit field and add degree
    fn pack_exponents(exps: [u8; 6]) -> u64 {
        let mut packed: u64 = 0;
        for (i, &e) in exps.iter().enumerate() {
            let shift = 40_u32.saturating_sub((8 * i) as u32);
            packed |= (e as u64) << shift;
        }
        let degree: u16 = exps.iter().map(|&x| x as u16).sum();
        packed | ((degree as u64) << 48)
    }

    // Helper: unpack the 6 exponents (ignores degree)
    fn unpack_exponents(packed: u64) -> [u8; 6] {
        let mut exps = [0u8; 6];
        for i in 0..6 {
            let shift = 40_u32.saturating_sub((8 * i) as u32);
            exps[i] = ((packed >> shift) & 0xFF) as u8;
        }
        exps
    } */

   /*  // Helper: per-field add of packed exponent vectors (returns packed with recomputed degree)
    fn add_packed(a: u64, b: u64) -> u64 {
        let ea = Term::unpack_exponents(a);
        let eb = Term::unpack_exponents(b);
        let mut er = [0u8; 6];
        for i in 0..6 {
            er[i] = ea[i].saturating_add(eb[i]);
        }
        Term::pack_exponents(er)
    }

    // Helper: per-field subtract (assumes b <= a component-wise)
    fn sub_packed(a: u64, b: u64) -> u64 {
        let ea = Term::unpack_exponents(a);
        let eb = Term::unpack_exponents(b);
        let mut er = [0u8; 6];
        for i in 0..6 {
            er[i] = ea[i].saturating_sub(eb[i]);
        }
        Term::pack_exponents(er)
    } */
    pub fn compare(&self, other: &Term) -> std::cmp::Ordering {
        let order = unsafe { TERM_ORDER };
        match order {
            TermOrder::Lex => {
                // Compare packed exponent fields directly (ignore degree)
                let a = self.exponents & 0x0000_FFFF_FFFF_FFFF;
                let b = other.exponents & 0x0000_FFFF_FFFF_FFFF;
                a.cmp(&b)
            }
            TermOrder::GrLex => {
                let self_degree = (self.exponents >> 48) & 0xFFFF;
                let other_degree = (other.exponents >> 48) & 0xFFFF;
                if self_degree != other_degree {
                    return self_degree.cmp(&other_degree);
                }
                // Compare packed exponents directly
                let a = self.exponents & 0x0000_FFFF_FFFF_FFFF;
                let b = other.exponents & 0x0000_FFFF_FFFF_FFFF;
                a.cmp(&b)
            }
            TermOrder::RevLex => {
                let self_degree = (self.exponents >> 48) & 0xFFFF;
                let other_degree = (other.exponents >> 48) & 0xFFFF;
                if self_degree != other_degree {
                    return self_degree.cmp(&other_degree);
                }
                // Reverse compare packed exponents
                let a = self.exponents & 0x0000_FFFF_FFFF_FFFF;
                let b = other.exponents & 0x0000_FFFF_FFFF_FFFF;
                b.cmp(&a) // Reverse order
            }
        }
    }

    pub fn can_reduce(&self, divisor: &Term) -> bool {
        let self_exponents = self.exponents & 0x0000_FFFF_FFFF_FFFF;
        let divisor_exponents = divisor.exponents & 0x0000_FFFF_FFFF_FFFF;

        for i in (0..48).step_by(8) {
            let self_exp = (self_exponents >> i) & 0xFF;
            let divisor_exp = (divisor_exponents >> i) & 0xFF;

            if divisor_exp > self_exp {
                return false;
            }
        }
        true
    }

    pub fn LCM(&self, other: &Term) -> u64 {
        let self_exponents = self.exponents & 0x0000_FFFF_FFFF_FFFF;
        let other_exponents = other.exponents & 0x0000_FFFF_FFFF_FFFF;

        let mut lcm_exponents: u64 = 0;
        let mut degree: u16 = 0;

        for i in (0..48).step_by(8) {
            let self_exp = (self_exponents >> i) & 0xFF;
            let other_exp = (other_exponents >> i) & 0xFF;

            let lcm_exp = self_exp.max(other_exp);

            lcm_exponents |= lcm_exp << i;
            degree += lcm_exp as u16; // Accumulate the degree
        }

        // Set the degree in the top 16 bits
        lcm_exponents |= (degree as u64) << 48;

        lcm_exponents
    }

    pub fn debug_print(&self) {
        print!("C:{} | ", self.coefficient);

        // Print the degree (top 16 bits) in hex
        let degree = (self.exponents >> 48) & 0xFFFF;
        print!("E(hex):{:04X} ", degree);

        // Print each 8-bit exponent in hex
        for i in (0..48).step_by(8).rev() {
            let exp = (self.exponents >> i) & 0xFF;
            print!("{:02X} ", exp);
        }
        println!(); // Newline after printing all exponents
    }
}


impl Eq for Term {}

impl Hash for Term {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.coefficient.hash(state);
        self.exponents.hash(state);
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct Polynomial {
    pub terms: Vec<Term>,
}

impl Eq for Polynomial {

}

impl Hash for Polynomial {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.terms.hash(state);
    }
}

fn mod_inverse(a: u64, m: u64) -> u64 {
    let (mut a, mut m) = (a as i64, m as i64);
    let (mut x0, mut x1) = (0i64, 1i64);
    let m0 = m;
    if m == 1 {
        return 0;
    }
    while a > 1 {
        let q = a / m;
        let t = m;
        m = a % m;
        a = t;
        let tmp = x0;
        x0 = x1 - q * x0;
        x1 = tmp;
    }
    if x1 < 0 {
        x1 += m0;
    }
    (x1 % m0) as u64
}

impl Polynomial {
    pub fn debug_print(&self) {
        println!("Polynomial:");
        for term in &self.terms {
            term.debug_print();
        }
    }

    pub fn new(mut terms: Vec<Term>) -> Self {
        // Sort terms by sort order
        terms.sort_by(|a, b| b.compare(a));
        terms.retain(|t| t.coefficient != 0); // Remove zero coefficient terms
       
        Polynomial { terms }
    }

    pub fn add(&self, other: &Polynomial) -> Polynomial {
        let modulus = unsafe { MODULUS };
        let mut result = self.terms.clone();
        for term in &other.terms {
            let mut found = false;
            for res_term in &mut result {
                if res_term.exponents == term.exponents {
                    res_term.coefficient = (res_term.coefficient + term.coefficient) % modulus;
                    found = true;
                    break;
                }
            }
            if !found {
                result.push(term.clone());
            }
        }
        // print results
        /*for term in &result {
            println!("Term: {:?}, Coefficient: {}", term.exponents, term.coefficient);
        }*/


        Polynomial::new(result)
    }

    pub fn subtract(&self, other: &Polynomial) -> Polynomial {
        let modulus = unsafe {MODULUS};        let mut result = self.terms.clone();
        for term in &other.terms {
            let mut found = false;
            for res_term in &mut result {
                if res_term.exponents == term.exponents {
                    res_term.coefficient = (modulus + res_term.coefficient - term.coefficient) % modulus;
                    found = true;
                    break;
                }
            }
            if !found {
                let mut neg_term = term.clone();
                neg_term.coefficient = (modulus + 0- term.coefficient) % modulus;
                result.push(neg_term);
            }
        }
        // print results
        /*for term in &result {
            println!("Term: {:?}, Coefficient: {}", term.exponents, term.coefficient);
        }*/

        Polynomial::new(result)
    }
    
    pub fn make_monic(&self) -> Polynomial {
        let modulus = unsafe {MODULUS};
        if self.terms.is_empty() { return self.clone(); }
        let lead_coeff = self.terms[0].coefficient;
        let inv = mod_inverse(lead_coeff, modulus);
        let new_terms = self.terms.iter().map(|t| Term {
            coefficient: (t.coefficient * inv) % modulus,
            exponents: t.exponents.clone(),
        }).collect();
        Polynomial::new(new_terms)
    }

   pub fn reduce(&self, divisors: &[Polynomial]) -> Polynomial {
        let modulus = unsafe {MODULUS};
        let mut result = self.clone(); // Start with the input polynomial
        let mut remainder: Vec<Term> = Vec::new();

        loop {
            let mut reduced = false;

            // Iterate over the divisors to reduce the leading term
            for divisor in divisors {
                // debug
                //println!("Using divisor:");
                //divisor.debug_print();
                if let Some(leading_term) = result.terms.first() {
                    if let Some(divisor_leading_term) = divisor.terms.first() {
                        // Check if the leading term can be reduced
                        if leading_term.can_reduce(&divisor_leading_term) {
                            let coefficient = leading_term.coefficient * mod_inverse(divisor_leading_term.coefficient, modulus) % modulus;
                            let exponents = leading_term.exponents - divisor_leading_term.exponents;

                            let reduction_term = Term {
                                coefficient,
                                exponents,
                            };

                            let scaled_divisor = divisor.multiply_by_term(&reduction_term);
                            result = result.subtract(&scaled_divisor);

                            reduced = true;
                            break; // Restart the loop after reducing
                        }
                    }
                }
            }

            // If no reduction was performed, break the loop
            if !reduced {
                if let Some(leading_term) = result.terms.first().cloned() {
                    //if(debug) { println!("No further reduction possible. Moving leading term {:?} to remainder.", leading_term); }
                    remainder.push(leading_term);
                    result.terms.remove(0);
                }
                else {
                    break;
                }
            }
        }

        result.terms.append(&mut remainder);
        Polynomial::new(result.terms)
    }

    pub fn multiply_by_term(&self, term: &Term) -> Polynomial {
        let modulus = unsafe {MODULUS};
        let terms = self
            .terms
            .iter()
            .map(|t| Term {
                coefficient: (t.coefficient * term.coefficient) % modulus,
                exponents: t.exponents + term.exponents,
            })
            .collect();

        Polynomial::new(terms)
    }

    pub fn s_polynomial(p1: &Polynomial, p2: &Polynomial) -> Polynomial {
        // Extract the leading terms of p1 and p2
        let leading_term_p1 = &p1.terms[0];
        let leading_term_p2 = &p2.terms[0];

        // Compute the LCM of the leading monomials' exponents
        let lcm_exponents = leading_term_p1.LCM(&leading_term_p2);

        let scale_factor_p1 = lcm_exponents - leading_term_p1.exponents;
        let scale_factor_p2 = lcm_exponents - leading_term_p2.exponents;

        let scaled_p1 = p1.multiply_by_term(&Term {
            coefficient: leading_term_p2.coefficient.clone(),
            exponents: scale_factor_p1,
        });

        let scaled_p2 = p2.multiply_by_term(&Term {
            coefficient: leading_term_p1.coefficient.clone(),
            exponents: scale_factor_p2,
        });

        scaled_p1.subtract(&scaled_p2)
    }

}

pub fn naive_grobner_basis(polynomials: Vec<Polynomial>) -> Vec<Polynomial> {
    let mut basis = polynomials.clone();
    let mut basis_set: HashSet<Polynomial> = HashSet::new();
    // print basis and polynomials
    for poly in &basis {
        //println!("{:?}", poly);
        basis_set.insert(poly.clone());
    }

    let mut pairs = Vec::<(usize, usize)>::new();
    for i in 0..basis.len() {
        for j in i + 1..basis.len() {
            pairs.push((i, j));
        }
    }
    //println!("Begin the experiment, {}", basis.len());
    while pairs.is_empty() == false {
        let (i, j) = pairs.remove(0);
        //println!("Processing pair ({}, {})", i, j);
        let s_poly = Polynomial::s_polynomial(&basis[i], &basis[j]);
        /*let mut debug = false;
        if(i == 0 && j == 7 || i == 3 && j == 4 || j ==4 && i==3) 
        { debug = true; println!("Debugging S-Polynomial for basis[{}] and basis[{}]", i, j); 
        // print basis
        for poly in &basis {
                println!("{:?}", poly);
            }
        }*/
        let reduced = s_poly.reduce(&basis);
        if !reduced.terms.is_empty() && !basis_set.contains(&reduced) {
            //println!("Adding new polynomial to basis."); 
            basis_set.insert(reduced.clone());
            let new_idx = basis.len();
            basis.push(reduced);
            pairs.extend((0..new_idx).map(|k| (k, new_idx)));
        }
    }
       

        //print basis with new lines separating each polynomial
        /*println!("New basis polynomials:");
        for poly in &basis {
            println!("{:?}", poly);
        }
        println!("End of iteration {}\n", i);*/
    

    //reduce basis by self
    let mut reduced_basis = Vec::new();
    for poly in &basis {
        // reduce poly by basis excluding itself
        let mut basis_excluding_self = basis.clone();
        basis_excluding_self.retain(|p| p != poly);
        let reduced = poly.reduce(&basis_excluding_self);
        if !reduced.terms.is_empty() && !reduced_basis.contains(&reduced) {
            reduced_basis.push(reduced.make_monic());
        }
    }
    reduced_basis
    //basis

    
}

pub fn are_bases_equivalent(set_a: Vec<Polynomial>, set_b: Vec<Polynomial>, modulus: u64) -> bool {
    // Check if all polynomials in set_a reduce to zero using set_b
    for poly in &set_a {
        let reduced = poly.reduce(&set_b);
        if !reduced.terms.is_empty() {
            return false; // Found a polynomial in set_a that does not reduce to zero
        }
    }

    // Check if all polynomials in set_b reduce to zero using set_a
    for poly in &set_b {
        let reduced = poly.reduce(&set_a);
        if !reduced.terms.is_empty() {
            return false; // Found a polynomial in set_b that does not reduce to zero
        }
    }

    // Both checks passed, the bases are equivalent
    true
}

fn main() {
    // let mode = 0 be for testing
    let mode = 0;

    if mode != 0 { 
        // arg1 = # of polynomials
        // arg2 = term order (0=Lex, 1=GrLex, 2=RevLex)
        let args: Vec<String> = std::env::args().collect();
        let mut rand = Lcg::new(12345, 1345, 16645, 1013904);
        let num_polynomials: usize = args.get(1).and_then(|s| s.parse().ok()).unwrap_or(3);
        let term_order: usize = args.get(2).and_then(|s| s.parse().ok()).unwrap_or(0);

        match term_order {
            0 => unsafe {TERM_ORDER = TermOrder::Lex},
            1 => unsafe {TERM_ORDER = TermOrder::GrLex},
            2 => unsafe {TERM_ORDER = TermOrder::RevLex},
            _ => unsafe {TERM_ORDER = TermOrder::Lex},
        }

        let mut input_basis = Vec::new();
        let modulus = 13 as u64;
        for _ in 0..num_polynomials {
            let mut terms = Vec::new();
            for _ in 0..3 { // always 3 terms per polynomial
                let coefficient = (rand.next_int() as u64 % modulus); // Coefficient between 0 and modulus-1
                let mut exponents = [0u8; 6];
                for i in 0..3 { // only working with 3 variables for now
                    exponents[i] = (rand.next_int() % 4) as u8; // Exponent between 0 and 3
                }
                terms.push(Term::from_exponents(coefficient, exponents));
            }
            input_basis.push(Polynomial::new(terms));
        }

        let basis = naive_grobner_basis(input_basis);
        println!("{}", basis.len());
        /*println!("Computed Grobner Basis Polynomials:");
        for poly in &basis {
            poly.debug_print();
            println!("---");
        }*/
        return;
    }
    else {
        let modulus = 7 as u64;
        let args: Vec<String> = std::env::args().collect();
        let n = args.get(1).and_then(|s| s.parse().ok()).unwrap_or(4);
        unsafe { TERM_ORDER = TermOrder::Lex; }

        if n == 4 {
            println!("Rust specialized finite coeff bitpacked cyclic 4");
            // Cyclic 4 system: x0+x1+x2+x3, x0x1+x1x2+x2x3+x3x0, x0x1x2+x1x2x3+x2x3x0+x3x0x1, x0x1x2x3-1
            // f1 = x0 + x1 + x2 + x3
            let p1 = Polynomial::new(vec![
                Term::from_exponents(1, [1, 0, 0, 0, 0, 0]),
                Term::from_exponents(1, [0, 1, 0, 0, 0, 0]),
                Term::from_exponents(1, [0, 0, 1, 0, 0, 0]),
                Term::from_exponents(1, [0, 0, 0, 1, 0, 0]),
            ]);

            // f2 = x0x1 + x1x2 + x2x3 + x3x0
            let p2 = Polynomial::new(vec![
                Term::from_exponents(1, [1, 1, 0, 0, 0, 0]),
                Term::from_exponents(1, [0, 1, 1, 0, 0, 0]),
                Term::from_exponents(1, [0, 0, 1, 1, 0, 0]),
                Term::from_exponents(1, [1, 0, 0, 1, 0, 0]),
            ]);

            // f3 = x0x1x2 + x1x2x3 + x2x3x0 + x3x0x1
            let p3 = Polynomial::new(vec![
                Term::from_exponents(1, [1, 1, 1, 0, 0, 0]),
                Term::from_exponents(1, [0, 1, 1, 1, 0, 0]),
                Term::from_exponents(1, [1, 0, 1, 1, 0, 0]),
                Term::from_exponents(1, [1, 1, 0, 1, 0, 0]),
            ]);

            // f4 = x0x1x2x3 - 1
            let p4 = Polynomial::new(vec![
                Term::from_exponents(1, [1, 1, 1, 1, 0, 0]),
                Term::from_exponents(modulus-1, [0, 0, 0, 0, 0, 0]),
            ]);
            let start = vec![p1,p2,p3,p4];
            for i in 0..10 {
                let basis = naive_grobner_basis(start.clone());
                println!("Iteration {}: complete", i);
                if i == 9 {
                    println!("Final Grobner Basis:");
                    for poly in &basis {
                        poly.debug_print();
                        println!("---");
                    }
                }

            }
        }
        else if n == 5 {
            println!("Rust specialized finite coeff bitpacked cyclic 5");
            // Cyclic 5
            let p1 = Polynomial::new(vec![
                Term::from_exponents(1, [1, 0, 0, 0, 0, 0]),
                Term::from_exponents(1, [0, 1, 0, 0, 0, 0]),
                Term::from_exponents(1, [0, 0, 1, 0, 0, 0]),
                Term::from_exponents(1, [0, 0, 0, 1, 0, 0]),
                Term::from_exponents(1, [0, 0, 0, 0, 1, 0]),
            ]);
            let p2 = Polynomial::new(vec![
                Term::from_exponents(1, [1, 1, 0, 0, 0, 0]),
                Term::from_exponents(1, [0, 1, 1, 0, 0, 0]),
                Term::from_exponents(1, [0, 0, 1, 1, 0, 0]),
                Term::from_exponents(1, [0, 0, 0, 1, 1, 0]),
                Term::from_exponents(1, [1, 0, 0, 0, 1, 0]),
            ]);
            let p3 = Polynomial::new(vec![
                Term::from_exponents(1, [1, 1, 1, 0, 0, 0]),
                Term::from_exponents(1, [0, 1, 1, 1, 0, 0]),
                Term::from_exponents(1, [0, 0, 1, 1, 1, 0]),
                Term::from_exponents(1, [1, 0, 0, 1, 1, 0]),
                Term::from_exponents(1, [1, 1, 0, 0, 1, 0]),
            ]);
            let p4 = Polynomial::new(vec![
                Term::from_exponents(1, [1, 1, 1, 1, 0, 0]),
                Term::from_exponents(1, [0, 1, 1, 1, 1, 0]),
                Term::from_exponents(1, [1, 0, 1, 1, 1, 0]),
                Term::from_exponents(1, [1, 1, 0, 1, 1, 0]),
                Term::from_exponents(1, [1, 1, 1, 0, 1, 0]),
            ]);
            let p5 = Polynomial::new(vec![
                Term::from_exponents(1, [1, 1, 1, 1, 1, 0]),
                Term::from_exponents(modulus-1, [0, 0, 0, 0, 0, 0]),
            ]);
            let start = vec![p1,p2,p3,p4,p5];
            for i in 0..10 {
                let basis = naive_grobner_basis(start.clone());
                println!("Iteration {}: complete", i);
                if i == 9 {
                    println!("Final Grobner Basis:");
                    for poly in &basis {
                        poly.debug_print();
                        println!("---");
                    }
                }
            }
        }
        else if n == 6 {
            println!("Rust specialized finite coeff bitpacked cyclic 6");
            // Cyclic 6
            let p1 = Polynomial::new(vec![
                Term::from_exponents(1, [1, 0, 0, 0, 0, 0]),
                Term::from_exponents(1, [0, 1, 0, 0, 0, 0]),
                Term::from_exponents(1, [0, 0, 1, 0, 0, 0]),
                Term::from_exponents(1, [0, 0, 0, 1, 0, 0]),
                Term::from_exponents(1, [0, 0, 0, 0, 1, 0]),
                Term::from_exponents(1, [0, 0, 0, 0, 0, 1]),
            ]);
             let p2 = Polynomial::new(vec![
                Term::from_exponents(1, [1, 1, 0, 0, 0, 0]),
                Term::from_exponents(1, [0, 1, 1, 0, 0, 0]),
                Term::from_exponents(1, [0, 0, 1, 1, 0, 0]),
                Term::from_exponents(1, [0, 0, 0, 1, 1, 0]),
                Term::from_exponents(1, [0, 0, 0, 0, 1, 1]),
                Term::from_exponents(1, [1 ,0 ,0 ,0 ,0 ,1]),
            ]);
             let p3 = Polynomial::new(vec![
                Term::from_exponents(1,[1 ,1 ,1 ,0 ,0 ,0]),
                Term::from_exponents(1,[0 ,1 ,1 ,1 ,0 ,0]),
                Term::from_exponents(1,[0 ,0 ,1 ,1 ,1 ,0]),
                Term::from_exponents(1,[0 ,0 ,0 ,1 ,1 ,1]),
                Term::from_exponents(1,[1 ,0 ,0 ,0 ,1 ,1]),
                Term::from_exponents(1,[1 ,1 ,0 ,0 ,0 ,1]),
            ]);
             let p4 = Polynomial::new(vec![
                Term::from_exponents(1,[1 ,1 ,1 ,1 ,0 ,0]),
                Term::from_exponents(1,[0 ,1 ,1 ,1 ,1 ,0]),
                Term::from_exponents(1,[0 ,0 ,1 ,1 ,1 ,1]),
                Term::from_exponents(1,[1 ,0 ,0 ,1 ,1 ,1]),
                Term::from_exponents(1,[1 ,1 ,0 ,0 ,1 ,1]),
                Term::from_exponents(1,[1 ,1 ,1 ,0 ,0 ,1]),
            ]);
             let p5 = Polynomial::new(vec![
                Term::from_exponents(1,[1 ,1 ,1 ,1 ,1 ,0]),
                Term::from_exponents(1,[0 ,1 ,1 ,1 ,1 ,1]),
                Term::from_exponents(1,[1 ,0 ,1 ,1 ,1 ,1]),
                Term::from_exponents(1,[1 ,1 ,0 ,1 ,1 ,1]),
                Term::from_exponents(1,[1 ,1 ,1 ,0 ,1 ,1]),
                Term::from_exponents(1,[1 ,1 ,1 ,1 ,0 ,1]),
            ]);
             let p6 = Polynomial::new(vec![
                Term::from_exponents(1,[1 ,1 ,1 ,1 ,1 ,1]),
                Term::from_exponents(modulus-1,[0 ,0 ,0 ,0 ,0 ,0]),
            ]);
            let start = vec![p1,p2,p3,p4,p5,p6];
            for i in 0..10 {
                let basis = naive_grobner_basis(start.clone());
                println!("Iteration {}: complete", i);
                if i == 9 {
                    println!("Final Grobner Basis:");
                    for poly in &basis {
                        poly.debug_print();
                        println!("---");
                    }
                }
            }
        }
        else if n == 7 {
            //bitpacked cant fit 7 variables
        }
    }
}
#[allow(dead_code)]
pub fn run_algorithm() {
    main();
}
