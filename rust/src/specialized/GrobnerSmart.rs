use std::collections::HashSet;

use std::hash::{Hash, Hasher};

use std::sync::OnceLock;

static TERM_ORDER: OnceLock<TermOrder> = OnceLock::new();

#[derive(Clone, Debug, PartialEq)]
pub struct Term {
    pub coefficient: f64,
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
    pub fn from_exponents(coefficient: f64, exponents: [u8; 6]) -> Self {
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
        let order = TERM_ORDER.get().expect("TERM_ORDER not initialized");
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
        println!("Coefficient: {}", self.coefficient);

        // Print the degree (top 16 bits) in hex
        let degree = (self.exponents >> 48) & 0xFFFF;
        print!("Degree: {:04X}, Exponents (hex): ", degree);

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
        self.coefficient.to_bits().hash(state);
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
        terms.retain(|t| t.coefficient != 0.0); // Remove zero coefficient terms
        // remove terms that are very close but not equal to 0 to handle floating point errors
        terms.retain(|t| (t.coefficient - 0.0).abs() > 1e-2);
        // round coefficients to 5 decimal places to handle floating point errors
        for term in &mut terms {
            term.coefficient = (term.coefficient * 1e5).round() / 1e5;
        }   
        Polynomial { terms }
    }

    pub fn add(&self, other: &Polynomial) -> Polynomial {
        let mut result = self.terms.clone();
        for term in &other.terms {
            let mut found = false;
            for res_term in &mut result {
                if res_term.exponents == term.exponents {
                    res_term.coefficient = res_term.coefficient + term.coefficient;
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
        let mut result = self.terms.clone();
        for term in &other.terms {
            let mut found = false;
            for res_term in &mut result {
                if res_term.exponents == term.exponents {
                    res_term.coefficient = res_term.coefficient - term.coefficient;
                    found = true;
                    break;
                }
            }
            if !found {
                let mut neg_term = term.clone();
                neg_term.coefficient = -term.coefficient;
                result.push(neg_term);
            }
        }
        // print results
        /*for term in &result {
            println!("Term: {:?}, Coefficient: {}", term.exponents, term.coefficient);
        }*/

        Polynomial::new(result)
    }
    

   pub fn reduce(&self, divisors: &[Polynomial]) -> Polynomial {
        let mut result = self.clone(); // Start with the input polynomial

        loop {
            let mut reduced = false;

            // Iterate over the divisors to reduce the leading term
            for divisor in divisors {
                if let Some(leading_term) = result.terms.first() {
                    if let Some(divisor_leading_term) = divisor.terms.first() {
                        // Check if the leading term can be reduced
                        if leading_term.can_reduce(&divisor_leading_term) {
                            let coefficient = leading_term.coefficient / divisor_leading_term.coefficient;
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
                break;
            }
        }

        Polynomial::new(result.terms)
    }

    pub fn multiply_by_term(&self, term: &Term) -> Polynomial {
        let terms = self
            .terms
            .iter()
            .map(|t| Term {
                coefficient: t.coefficient * term.coefficient,
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
            coefficient: 1.0,
            exponents: scale_factor_p1,
        });

        let scaled_p2 = p2.multiply_by_term(&Term {
            coefficient: 1.0,
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
        poly.debug_print();
        println!("---");
    }

    println!("Begin the experiment, {}", basis.len());
    for i in 0..10 { // This is *supposed* to go until no new polynomials are added, but for now just do 3 iterations
        let basis_len = basis.len();
        let mut added = false;
        for i in 0..basis_len {
            for j in i + 1..basis_len {
                let s_poly = Polynomial::s_polynomial(&basis[i], &basis[j]);
                let reduced = s_poly.reduce(&basis);
                //print basis[i], basis[j], s_poly, reduced
                //println!("Basis 1: {:?} | Basis 2: {:?} | S-Polynomial: {:?}", basis[i], basis[j], s_poly);
                //println!("Reduced: {:?}", reduced);
                if !reduced.terms.is_empty() && !basis_set.contains(&reduced) {
                    //println!("Adding new polynomial to basis.");
                    basis_set.insert(reduced.clone());
                    basis.push(reduced);
                    added = true;
                }
                else {
                    //println!("Reduced polynomial is zero or already in basis, skipping.");
                }
            }
        }

        if !added {
            break;
        }

        //print basis with new lines separating each polynomial
        //println!("New basis polynomials:");
        for poly in &basis {
            //println!("{:?}", poly);
        }
        println!("End of iteration {}\n", i);
    }

    //reduce basis by self
    let mut reduced_basis = Vec::new();
    for poly in &basis {
        // reduce poly by basis excluding itself
        let mut basis_excluding_self = basis.clone();
        basis_excluding_self.retain(|p| p != poly);
        let reduced = poly.reduce(&basis_excluding_self);
        if !reduced.terms.is_empty() && !reduced_basis.contains(&reduced) {
            reduced_basis.push(reduced);
        }
    }

    /* I don't think this is necessary, but it does make for a nicer output
    //if leading term of polynomial has negative coefficient flip all signs
    for poly in &mut reduced_basis {
        if let Some(leading_term) = poly.terms.first() {
            if leading_term.coefficient < 0.0 {
                for term in &mut poly.terms {
                    term.coefficient = -term.coefficient;
                }
            }
        }
    }*/

    reduced_basis
}

pub fn are_bases_equivalent(set_a: Vec<Polynomial>, set_b: Vec<Polynomial>) -> bool {
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
    // 1 for s_polynomial, 2 for add, 3 for subtract, 4 for reduce, 5 for testing hashes, else grobner basis
    let test = 0;
     //Lex, GrLex, RevLex
    TERM_ORDER.set(TermOrder::Lex).expect("TERM_ORDER already initialized");

    // x*y - z
    let p1 = Polynomial::new(vec![
        Term::from_exponents(1.0, [1, 1, 0, 0, 0, 0]),
        Term::from_exponents(-1.0, [0, 0, 1, 0, 0, 0]),
    ]);
    // x^2 + y^2 -1
    let p2 = Polynomial::new(vec![
        Term::from_exponents(1.0, [2, 0, 0, 0, 0, 0]),
        Term::from_exponents(1.0, [0, 2, 0, 0, 0, 0]),
        Term::from_exponents(-1.0, [0, 0, 0, 0, 0, 0]),
    ]);
      
    // x+y+z
    let p3 = Polynomial::new(vec![
        Term::from_exponents(1.0, [1, 0, 0, 0, 0, 0]),
        Term::from_exponents(1.0, [0, 1, 0, 0, 0, 0]),
        Term::from_exponents(1.0, [0, 0, 1, 0, 0, 0]),
    ]);

    // x + y + z
    let p4 = Polynomial::new(vec![
        Term::from_exponents(1.0, [1, 0, 0, 0, 0, 0]),
        Term::from_exponents(1.0, [0, 1, 0, 0, 0, 0]),
        Term::from_exponents(1.0, [0, 0, 1, 0, 0, 0]),
    ]);

    println!("Begin the experiment");


 
    let basis = naive_grobner_basis(vec![p1, p2]);
    // copy basis
    let mut copied_basis = basis.clone();
    println!("Final Grobner Basis:");
    for poly in basis {
        poly.debug_print();
        println!("---");
    }
    

    //  x^2 + y^2 -1
    let test_poly = Polynomial::new(vec![
        Term::from_exponents(1.0, [2, 0, 0, 0, 0, 0]),
        Term::from_exponents(1.0, [0, 2, 0, 0, 0, 0]),
        Term::from_exponents(-1.0, [0, 0, 0, 0, 0, 0]),
    ]);

    // x*y -z
    let test_poly_2 = Polynomial::new(vec![
        Term::from_exponents(1.0, [1, 1, 0, 0, 0, 0]),
        Term::from_exponents(-1.0, [0, 0, 1, 0, 0, 0]),
        
    ]);

    // x*z +y^3 - y
    let test_poly_3 = Polynomial::new(vec![
        Term::from_exponents(1.0, [1, 0, 1, 0, 0, 0]),
        Term::from_exponents(1.0, [0, 3, 0, 0, 0, 0]),
        Term::from_exponents(-1.0, [0, 1, 0, 0, 0, 0]),
        
    ]);

    // y^4 - y^2 + z^2
    let test_poly_4 = Polynomial::new(vec![
        Term::from_exponents(1.0, [0, 4, 0, 0, 0, 0]),
        Term::from_exponents(-1.0, [0, 2, 0, 0, 0, 0]),
        Term::from_exponents(1.0, [0, 0, 2, 0, 0, 0]),
    ]);

    let test_basis = vec![test_poly, test_poly_2, test_poly_3, test_poly_4];

    // print test_basis
    println!("Test Basis Polynomials:");
    for poly in &test_basis {
        poly.debug_print();
        println!("---");
    }
    let is_equivalent = are_bases_equivalent(copied_basis, test_basis);
    println!("Are the computed basis and test basis equivalent? {}", is_equivalent);
    
}