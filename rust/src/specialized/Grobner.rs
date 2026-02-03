use std::collections::HashSet;

use std::hash::{Hash, Hasher};

use std::sync::OnceLock;

static TERM_ORDER: OnceLock<TermOrder> = OnceLock::new();

use rust::helpers::lcg::Lcg;

#[derive(Clone, Debug, PartialEq)]
pub struct Term {
    pub coefficient: f64,
    pub exponents: Vec<usize>, // Exponents for each variable
}

#[derive(Copy, Clone, Debug, PartialEq)]
pub enum TermOrder {
    Lex,
    GrLex,
    RevLex
}

impl Term {
    pub fn compare(&self, other: &Term) -> std::cmp::Ordering {
        let order = TERM_ORDER.get().expect("TERM_ORDER not initialized");
        match order {
            TermOrder::Lex => self.exponents.cmp(&other.exponents), // Lexicographic order
            TermOrder::GrLex => { // Graded lexicographic order
                let self_degree: usize = self.exponents.iter().sum();
                let other_degree: usize = other.exponents.iter().sum();
                self_degree.cmp(&other_degree).then_with(|| self.exponents.cmp(&other.exponents))
            }
            TermOrder::RevLex => { // Reverse lexicographic order
                let self_degree: usize = self.exponents.iter().sum();
                let other_degree: usize = other.exponents.iter().sum();
                self_degree
                    .cmp(&other_degree)
                    .then_with(|| self.exponents.iter().rev().cmp(other.exponents.iter().rev()))
            }
        }
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
    pub fn new(mut terms: Vec<Term>) -> Self {
        // Sort terms by sort order
        terms.sort_by(|a, b| b.compare(a));
        //terms.retain(|t| t.coefficient != 0.0); // Remove zero coefficient terms
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
                        if leading_term.exponents.iter().zip(&divisor_leading_term.exponents).all(|(a, b)| a >= b) {
                            // print leading_term and divisor_leading_term
                            //println!("Leading Term: {:?}, Divisor Leading Term: {:?}", leading_term, divisor_leading_term);

                            // Compute the reduction factor
                            let coefficient = leading_term.coefficient / divisor_leading_term.coefficient;
                            let exponents: Vec<usize> = leading_term
                                .exponents
                                .iter()
                                .zip(&divisor_leading_term.exponents)
                                .map(|(a, b)| a - b)
                                .collect();

                            // Create the reduction term
                            let reduction_term = Term {
                                coefficient,
                                exponents,
                            };

                            // Subtract the scaled divisor from the result
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
                exponents: t
                    .exponents
                    .iter()
                    .zip(&term.exponents)
                    .map(|(a, b)| a + b)
                    .collect(),
            })
            .collect();

        Polynomial::new(terms)
    }

    pub fn s_polynomial(p1: &Polynomial, p2: &Polynomial) -> Polynomial {
        // Compute the LCM of the leading monomials' exponents
        let mut lcm_exponents = vec![0; p1.terms[0].exponents.len()];
        for i in 0..lcm_exponents.len() {
            lcm_exponents[i] = p1.terms[0].exponents[i].max(p2.terms[0].exponents[i]);
        }
        
        // Scale p1 to the LCM
        let scale_factor_p1 = lcm_exponents
            .iter()
            .zip(&p1.terms[0].exponents)
            .map(|(lcm, exp)| lcm - exp)
            .collect::<Vec<_>>();

        let scaled_p1 = Polynomial::new(
            p1.terms
                .iter()
                .map(|term| Term {
                    coefficient: term.coefficient,
                    exponents: term
                        .exponents
                        .iter()
                        .zip(&scale_factor_p1)
                        .map(|(exp, scale)| exp + scale)
                        .collect(),
                })
                .collect(),
        );

        // Scale p2 to the LCM
        let scale_factor_p2 = lcm_exponents
            .iter()
            .zip(&p2.terms[0].exponents)
            .map(|(lcm, exp)| lcm - exp)
            .collect::<Vec<_>>();

        let scaled_p2 = Polynomial::new(
            p2.terms
                .iter()
                .map(|term| Term {
                    coefficient: term.coefficient,
                    exponents: term
                        .exponents
                        .iter()
                        .zip(&scale_factor_p2)
                        .map(|(exp, scale)| exp + scale)
                        .collect(),
                })
                .collect(),
        );

        // Subtract the scaled polynomials
        scaled_p1.subtract(&scaled_p2)
    }

}

pub fn naive_grobner_basis(polynomials: Vec<Polynomial>) -> Vec<Polynomial> {
    let mut basis = polynomials.clone();
    let mut basis_set: HashSet<Polynomial> = HashSet::new();
    // print basis and polynomials
    /*for poly in &basis {
        println!("{:?}", poly);
        basis_set.insert(poly.clone());
    }*/

    //println!("Begin the experiment, {}", basis.len());
    loop { 
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
                    //println!("Basis 1: {:?} | Basis 2: {:?} | S-Polynomial: {:?}", basis[i], basis[j], s_poly);
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
        /*println!("New basis polynomials:");
        for poly in &basis {
            println!("{:?}", poly);
        }*/
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
    // let mode = 0 be for testing
    let mode = 0;
    if mode != 0 {
        // arg1 = # of polynomials
        // arg2 = term order (0=Lex, 1=GrLex, 2=RevLex)
        let args: Vec<String> = std::env::args().collect();
        let mut rand = Lcg::new(12345, 1345, 65, 17);
        let num_polynomials: usize = args.get(1).and_then(|s| s.parse().ok()).unwrap_or(3);
        let term_order: usize = args.get(2).and_then(|s| s.parse().ok()).unwrap_or(0);

        match term_order {
            0 => TERM_ORDER.set(TermOrder::Lex).expect("TERM_ORDER already initialized"),
            1 => TERM_ORDER.set(TermOrder::GrLex).expect("TERM_ORDER already initialized"),
            2 => TERM_ORDER.set(TermOrder::RevLex).expect("TERM_ORDER already initialized"),
            _ => TERM_ORDER.set(TermOrder::Lex).expect("TERM_ORDER already initialized"),
        }

        let mut input_basis = Vec::new();
        for _ in 0..num_polynomials {
            let mut terms = Vec::new();
            for _ in 0..3 { // always 3 terms per polynomial
                let coefficient = (rand.next_double() * 2.0 - 1.0);
                // only working with 3 variables for now
                let exponents = vec![(rand.next_int() % 4) as usize, (rand.next_int() % 4) as usize, (rand.next_int() % 4) as usize];
                terms.push(Term {
                    coefficient,
                    exponents,
                });
            }
            input_basis.push(Polynomial::new(terms));
        }

        let basis = naive_grobner_basis(input_basis);
        println!("{}", basis.len());
        /*println!("Computed Grobner Basis Polynomials:");
        for poly in &basis {
            println!("{:?}", poly);
        }*/
        return;
    }
    else {    
        // Cyclic-4 benchmark
        TERM_ORDER.set(TermOrder::Lex).expect("TERM_ORDER already initialized");

        // f1 = x0 + x1 + x2 + x3
        let p1 = Polynomial::new(vec![
            Term { coefficient: 1.0, exponents: vec![1, 0, 0, 0] },
            Term { coefficient: 1.0, exponents: vec![0, 1, 0, 0] },
            Term { coefficient: 1.0, exponents: vec![0, 0, 1, 0] },
            Term { coefficient: 1.0, exponents: vec![0, 0, 0, 1] },
        ]);

        // f2 = x0*x1 + x1*x2 + x2*x3 + x3*x0
        let p2 = Polynomial::new(vec![
            Term { coefficient: 1.0, exponents: vec![1, 1, 0, 0] },
            Term { coefficient: 1.0, exponents: vec![0, 1, 1, 0] },
            Term { coefficient: 1.0, exponents: vec![0, 0, 1, 1] },
            Term { coefficient: 1.0, exponents: vec![1, 0, 0, 1] },
        ]);

        // f3 = x0*x1*x2 + x1*x2*x3 + x2*x3*x0 + x3*x0*x1
        let p3 = Polynomial::new(vec![
            Term { coefficient: 1.0, exponents: vec![1, 1, 1, 0] },
            Term { coefficient: 1.0, exponents: vec![0, 1, 1, 1] },
            Term { coefficient: 1.0, exponents: vec![1, 0, 1, 1] },
            Term { coefficient: 1.0, exponents: vec![1, 1, 0, 1] },
        ]);

        // f4 = x0*x1*x2*x3 - 1
        let p4 = Polynomial::new(vec![
            Term { coefficient: 1.0, exponents: vec![1, 1, 1, 1] },
            Term { coefficient: -1.0, exponents: vec![0, 0, 0, 0] },
        ]);

        println!("Computing Grobner basis for Cyclic-4...");

        let basis = naive_grobner_basis(vec![p1, p2, p3, p4]);
        
        println!("Final Grobner Basis for Cyclic-4:");
        println!("Number of polynomials in basis: {}", basis.len());
        for (i, poly) in basis.iter().enumerate() {
            println!("Polynomial {}: {:?}", i + 1, poly);
        }

        // q1 = x0 + x1 + x2 + x3
        let q1 = Polynomial::new(vec![
            Term { coefficient: 1.0, exponents: vec![1, 0, 0, 0] },
            Term { coefficient: 1.0, exponents: vec![0, 1, 0, 0] },
            Term { coefficient: 1.0, exponents: vec![0, 0, 1, 0] },
            Term { coefficient: 1.0, exponents: vec![0, 0, 0, 1] },
        ]);

        // q2 = x1^2 + 2x1x3 + x3^2
        let q2 = Polynomial::new(vec![
            Term { coefficient: 1.0, exponents: vec![0, 2, 0, 0] },
            Term { coefficient: 2.0, exponents: vec![0, 1, 0, 1] },
            Term { coefficient: 1.0, exponents: vec![0, 0, 2, 0] },
        ]);

        // q3 = x1x2 - x1x3 + x2^2*x3^4 +x2x3 - 2x3^2
        let q3 = Polynomial::new(vec![
            Term { coefficient: 1.0, exponents: vec![0, 1, 1, 0] },
            Term { coefficient: -1.0, exponents: vec![0, 1, 0, 1] },
            Term { coefficient: 1.0, exponents: vec![0, 0, 2, 4] },
            Term { coefficient: 1.0, exponents: vec![0, 0, 1, 1] },
            Term { coefficient: -2.0, exponents: vec![0, 0, 0, 2] },
        ]);

        // q4 = x1*x3^4 -x1 + x3^5 - x3
        let q4 = Polynomial::new(vec![
            Term { coefficient: 1.0, exponents: vec![0, 1, 0, 4] },
            Term { coefficient: -1.0, exponents: vec![0, 1, 0, 0] },
            Term { coefficient: 1.0, exponents: vec![0, 0, 0, 5] },
            Term { coefficient: -1.0, exponents: vec![0, 0, 0, 1] },
        ]);

        // q5 = x2^3*x3^2 +x2^2*x3^3 - x2 - x3
        let q5 = Polynomial::new(vec![
            Term { coefficient: 1.0, exponents: vec![0, 0, 3, 2] },
            Term { coefficient: 1.0, exponents: vec![0, 0, 2, 3] },
            Term { coefficient: -1.0, exponents: vec![0, 0, 1, 0] },
            Term { coefficient: -1.0, exponents: vec![0, 0, 0, 1] },
        ]);

        // q6 = x2^2*x3^6 - x2^2*x3^2 -x3^4 + 1
        let q6 = Polynomial::new(vec![
            Term { coefficient: 1.0, exponents: vec![0, 0, 2, 6] },
            Term { coefficient: -1.0, exponents: vec![0, 0, 2, 2] },
            Term { coefficient: -1.0, exponents: vec![0, 0, 0, 4] },
            Term { coefficient: 1.0, exponents: vec![0, 0, 0, 0] },
        ]);

        let benchmark_basis = vec![q1, q2, q3, q4, q5, q6];
        let equivalent = are_bases_equivalent(basis, benchmark_basis);
        if equivalent {
            println!("The computed Grobner basis is equivalent to the benchmark basis.");
        } else {
            println!("The computed Grobner basis is NOT equivalent to the benchmark basis.");   
        }

    }
}
