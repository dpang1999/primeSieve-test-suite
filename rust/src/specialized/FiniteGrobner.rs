use std::collections::HashSet;

use std::hash::{Hash, Hasher};

static mut TERM_ORDER: TermOrder = TermOrder::Lex; // default to lex order, can be set to GrLex or RevLex as well

use rust::helpers::lcg::Lcg;

#[derive(Clone, Debug, PartialEq)]
pub struct Term {
    pub coefficient: u32,
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
        let order = unsafe { TERM_ORDER };
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

fn mod_inverse(a: u32, m: u32) -> u32 {
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
    (x1 % m0) as u32
}
impl Polynomial {
    pub fn new(mut terms: Vec<Term>) -> Self {
        terms.sort_by(|a, b| b.compare(a));
        terms.retain(|t| t.coefficient != 0); // Remove zero coefficient terms
        // remove terms that are very close but not equal to 0 to handle floating point errors
        //terms.retain(|t| (t.coefficient - 0.0).abs() > 1e-2);
        // round coefficients to 5 decimal places to handle floating point errors
        /*for term in &mut terms {
            term.coefficient = (term.coefficient * 1e5).round() / 1e5;
        } */  
        Polynomial { terms }
    }
    

    pub fn add(&self, other: &Polynomial, modulus: u32) -> Polynomial {
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

    pub fn subtract(&self, other: &Polynomial, modulus: u32) -> Polynomial {
        let mut result = self.terms.clone();
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
                neg_term.coefficient = (modulus + 0 - term.coefficient) % modulus;
                result.push(neg_term);
            }
        }
        // print results
        /*for term in &result {
            println!("Term: {:?}, Coefficient: {}", term.exponents, term.coefficient);
        }*/

        Polynomial::new(result)
    }

    pub fn make_monic(&self, modulus: u32) -> Polynomial {
        if self.terms.is_empty() { return self.clone(); }
        let lead_coeff = self.terms[0].coefficient;
        let inv = mod_inverse(lead_coeff, modulus);
        let new_terms = self.terms.iter().map(|t| Term {
            coefficient: (t.coefficient * inv) % modulus,
            exponents: t.exponents.clone(),
        }).collect();
        Polynomial::new(new_terms)
    }
    

   pub fn reduce(&self, divisors: &[Polynomial], modulus: u32) -> Polynomial {
        let mut result = self.clone(); // Start with the input polynomial
        let mut remainder: Vec<Term> = Vec::new();

        loop {
            let mut reduced = false;

            // Iterate over the divisors to reduce the leading term
            for divisor in divisors {
                if let Some(leading_term) = result.terms.first() {
                    if let Some(divisor_leading_term) = divisor.terms.first() {
                        // Check if the leading term can be reduced
                        if leading_term.exponents.iter().zip(&divisor_leading_term.exponents).all(|(a, b)| a >= b) {
                            // print leading_term and divisor_leading_term
                            //if (debug) {println!("Current f: {:?}\nLeading Term: {:?}, Divisor Leading Term: {:?}", result,leading_term, divisor_leading_term);}

                            // Compute the reduction factor
                            //let coefficient = leading_term.coefficient / divisor_leading_term.coefficient;
                            let coefficient = (leading_term.coefficient as u64 * mod_inverse(divisor_leading_term.coefficient, modulus) as u64 % modulus as u64) as u32;
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

                            //if (debug) {println!("Reduction Term: {:?}, Coefficient: {} * Dibisor: {:?}", reduction_term.exponents, reduction_term.coefficient, divisor);}
                            //println!("Divisor: {:?}", divisor);

                            // Subtract the scaled divisor from the result
                            let scaled_divisor = divisor.multiply_by_term(&reduction_term, modulus);
                            result = result.subtract(&scaled_divisor, modulus);

                            //if (debug) {println!("Result after reduction: {:?}\n", result);}

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

    pub fn multiply_by_term(&self, term: &Term, modulus: u32) -> Polynomial {
        let terms = self
            .terms
            .iter()
            .map(|t| Term {
                coefficient: (t.coefficient * term.coefficient) % modulus,
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

    pub fn s_polynomial(p1: &Polynomial, p2: &Polynomial, modulus: u32) -> Polynomial {
        // Compute the LCM of the leading monomials' exponents
        let mut lcm_exponents = vec![0; p1.terms[0].exponents.len()];
        for i in 0..lcm_exponents.len() {
            lcm_exponents[i] = p1.terms[0].exponents[i].max(p2.terms[0].exponents[i]);
        }

        // Compute exponent shifts for each polynomial
        let shift_p1 = lcm_exponents.iter().zip(&p1.terms[0].exponents).map(|(lcm, exp)| lcm - exp).collect::<Vec<_>>();
        let shift_p2 = lcm_exponents.iter().zip(&p2.terms[0].exponents).map(|(lcm, exp)| lcm - exp).collect::<Vec<_>>();

        // Compute scaling factors for coefficients (modular inverses)
        let a = p1.terms[0].coefficient;
        let b = p2.terms[0].coefficient;
        //let b_inv = mod_inverse(b, modulus);
        //let a_inv = mod_inverse(a, modulus);

        // S = (b * x^shift_p1 * p1 - a * x^shift_p2 * p2) mod modulus
        let scaled_p1 = Polynomial::new(
            p1.terms.iter().map(|term| Term {
                coefficient: (b * term.coefficient) % modulus,
                exponents: term.exponents.iter().zip(&shift_p1).map(|(exp, shift)| exp + shift).collect(),
            }).collect()
        );
        let scaled_p2 = Polynomial::new(
            p2.terms.iter().map(|term| Term {
                coefficient: (a * term.coefficient) % modulus,
                exponents: term.exponents.iter().zip(&shift_p2).map(|(exp, shift)| exp + shift).collect(),
            }).collect()
        );
        scaled_p1.subtract(&scaled_p2, modulus)
    }

}

pub fn naive_grobner_basis(polynomials: Vec<Polynomial>, modulus: u32) -> Vec<Polynomial> {
    let mut basis = polynomials.clone();
    let mut basis_set: HashSet<Polynomial> = HashSet::new();
    // print basis and polynomials
    for poly in &basis {
        //println!("{:?}", poly);
        basis_set.insert(poly.clone());
    }

    let mut processed_pairs = HashSet::<(usize, usize)>::new();
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
        processed_pairs.insert((i,j));
        let s_poly = Polynomial::s_polynomial(&basis[i], &basis[j], modulus);
        /*let mut debug = false;
        if(i == 0 && j == 7 || i == 3 && j == 4 || j ==4 && i==3) 
        { debug = true; println!("Debugging S-Polynomial for basis[{}] and basis[{}]", i, j); 
        // print basis
        for poly in &basis {
                println!("{:?}", poly);
            }
        }*/
        let reduced = s_poly.reduce(&basis, modulus);
        if !reduced.terms.is_empty() && !basis_set.contains(&reduced) {
            //println!("Adding new polynomial to basis."); 
            basis_set.insert(reduced.clone());
      
            basis.push(reduced);
            let new_idx = basis.len() - 1;
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
        let reduced = poly.reduce(&basis_excluding_self, modulus);
        if !reduced.terms.is_empty() && !reduced_basis.contains(&reduced) {
            reduced_basis.push(reduced.make_monic(modulus));
        }
    }
    reduced_basis
    //basis

    
}

pub fn are_bases_equivalent(set_a: Vec<Polynomial>, set_b: Vec<Polynomial>, modulus: u32) -> bool {
    // Check if all polynomials in set_a reduce to zero using set_b
    let mut all_ok = true;
    for (i, poly) in set_a.iter().enumerate() {
        let reduced = poly.reduce(&set_b, modulus);
        if !reduced.terms.is_empty() {
            println!("[DEBUG] set_a[{}] does not reduce to zero with set_b:", i);
            println!("  Original: {:?}", poly);
            println!("  Reduced:  {:?}", reduced);
            all_ok = false;
        }
    }

    // Check if all polynomials in set_b reduce to zero using set_a
    for (i, poly) in set_b.iter().enumerate() {
        let reduced = poly.reduce(&set_a, modulus);
        if !reduced.terms.is_empty() {
            println!("[DEBUG] set_b[{}] does not reduce to zero with set_a:", i);
            println!("  Original: {:?}", poly);
            println!("  Reduced:  {:?}", reduced);
            all_ok = false;
        }
    }

    if !all_ok {
        println!("[DEBUG] Bases are NOT equivalent.");
    } else {
        println!("[DEBUG] Bases are equivalent.");
    }
    all_ok
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
            0 => unsafe { TERM_ORDER = TermOrder::Lex; },
            1 => unsafe { TERM_ORDER = TermOrder::GrLex; },
            2 => unsafe { TERM_ORDER = TermOrder::RevLex; },
            _ => unsafe { TERM_ORDER = TermOrder::Lex; },
        }

        let mut input_basis = Vec::new();
        let modulus = 13;
        for _ in 0..num_polynomials {
            let mut terms = Vec::new();
            for _ in 0..3 { // always 3 terms per polynomial
                let coefficient = (rand.next_int() % (modulus as i32)) as u32;
                // only working with 3 variables for now
                let exponents = vec![(rand.next_int() % 4) as usize, (rand.next_int() % 4) as usize, (rand.next_int() % 4) as usize];
                terms.push(Term {
                    coefficient,
                    exponents,
                });
            }
            input_basis.push(Polynomial::new(terms));
        }

        let basis = naive_grobner_basis(input_basis, modulus);
        println!("{}", basis.len());
        println!("Computed Grobner Basis Polynomials:");
        for poly in &basis {
            println!("{:?}", poly);
        }
        return;
    }
    else {
        // 1 for s_polynomial, 2 for add, 3 for subtract, 4 for reduce, 5 for testing hashes, else grobner basis
        let test = 0;
        let modulus = 7;
        //Lex, GrLex, RevLex
        unsafe { TERM_ORDER = TermOrder::Lex; }

        // f1 = x0 + x1 + x2 + x3 + x4
        let p1 = Polynomial::new(vec![
            Term { coefficient: 1, exponents: vec![1, 0, 0, 0, 0] },
            Term { coefficient: 1, exponents: vec![0, 1, 0, 0, 0] },
            Term { coefficient: 1, exponents: vec![0, 0, 1, 0, 0] },
            Term { coefficient: 1, exponents: vec![0, 0, 0, 1, 0] },
            Term { coefficient: 1, exponents: vec![0, 0, 0, 0, 1] },
        ]);

        // f2 = x0x1 + x1x2 + x2x3 + x3x4 + x4x0
        let p2 = Polynomial::new(vec![
            Term { coefficient: 1, exponents: vec![1, 1, 0, 0, 0] },
            Term { coefficient: 1, exponents: vec![0, 1, 1, 0, 0] },
            Term { coefficient: 1, exponents: vec![0, 0, 1, 1, 0] },
            Term { coefficient: 1, exponents: vec![0, 0, 0, 1, 1] },
            Term { coefficient: 1, exponents: vec![1, 0, 0, 0, 1] },
        ]);

        // f3 = x0x1x2 + x1x2x3 + x2x3x4 + x3x4x0 + x4x0x1
        let p3 = Polynomial::new(vec![
            Term { coefficient: 1, exponents: vec![1, 1, 1, 0, 0] },
            Term { coefficient: 1, exponents: vec![0, 1, 1, 1, 0] },
            Term { coefficient: 1, exponents: vec![0, 0, 1, 1, 1] },
            Term { coefficient: 1, exponents: vec![1, 0, 0, 1, 1] },
            Term { coefficient: 1, exponents: vec![1, 1, 0, 0, 1] },
        ]);

        // f4 = x0x1x2x3 + x1x2x3x4 + x2x3x4x0 + x3x4x0x1 + x4x0x1x2 - 1
        let p4 = Polynomial::new(vec![
            Term { coefficient: 1, exponents: vec![1, 1, 1, 1, 0] },
            Term { coefficient: 1, exponents: vec![0, 1, 1, 1, 1] },
            Term { coefficient: 1, exponents: vec![1, 0, 1, 1, 1] },
            Term { coefficient: 1, exponents: vec![1, 1, 0, 1, 1] },
            Term { coefficient: 1, exponents: vec![1, 1, 1, 0, 1] },
        ]);

        // f4 = x0*x1*x2*x3*x4 -1
        let p5 = Polynomial::new(vec![
            Term { coefficient: 1, exponents: vec![1, 1, 1, 1, 1] },
            Term { coefficient: modulus-1, exponents: vec![0, 0, 0, 0, 0] },
        ]);

        println!("Begin the experiment");


        if test == 1 {
            let s_poly = Polynomial::s_polynomial(&p1, &p2, modulus);

            let expected = Polynomial::new(vec![
                Term {
                    coefficient: modulus - 1,    
                    exponents: vec![0, 2], // Term: -y^2
                },
                Term {
                    coefficient: 1,
                    exponents: vec![1, 0], // Term: x
                },
            ]);

            assert_eq!(s_poly, expected);
        } else if test == 2 {
            let sum = p1.add(&p2, modulus);
            println!("Sum: {:?}", sum);

            let expected = Polynomial::new(vec![
                Term {
                    coefficient: 1,
                    exponents: vec![2, 0],
                },
                Term {
                    coefficient: 1,
                    exponents: vec![1, 1],
                },
                Term {
                    coefficient: modulus - 1,
                    exponents: vec![0, 1],
                },
                Term {
                    coefficient: modulus - 1,
                    exponents: vec![0, 0],
                },
            ]);

            assert_eq!(sum, expected);
        } else if test == 3 {
            let difference = p1.subtract(&p2, modulus);
            println!("Difference: {:?}", difference);

            let expected = Polynomial::new(vec![
                Term {
                    coefficient: 1,
                    exponents: vec![2, 0],
                },
                Term {
                    coefficient: modulus - 1,
                    exponents: vec![1, 1],
                },
                Term {
                    coefficient: modulus - 1,
                    exponents: vec![0, 1],
                },
                Term {
                    coefficient: 1,
                    exponents: vec![0, 0],
                },
            ]);

            assert_eq!(difference, expected);
        }
        else if test == 4 {
            // x^3 + y^3 + z^3
            let p1 = Polynomial::new(vec![
                Term {
                    coefficient: 1,
                    exponents: vec![3, 0, 0], // x^3
                },
                Term {
                    coefficient: 1,
                    exponents: vec![0, 3, 0], // y^3
                },
                Term {
                    coefficient: 1,
                    exponents: vec![0, 0, 3], // z^3
                },
            ]);
            // xy + yz + xz
            let p2 = Polynomial::new(vec![
                Term {
                    coefficient: 1,
                    exponents: vec![1, 1, 0], // xy
                },
                Term {
                    coefficient: 1,
                    exponents: vec![0, 1, 1], // yz
                },
                Term {
                    coefficient: 1,
                    exponents: vec![1, 0, 1], // xz
                },
            ]);
            // x+y+z
            let p3 = Polynomial::new(vec![
                Term {
                    coefficient: 1,
                    exponents: vec![1, 0, 0], // x
                },
                Term {
                    coefficient: 1,
                    exponents: vec![0, 1, 0], // y
                },
                Term {
                    coefficient: 1,
                    exponents: vec![0, 0, 1], // z
                },
            ]);
            let basis = vec![p1.clone(), p2.clone(), p3.clone()];

            // xy^2 - y
            let divisor = Polynomial::new(vec![
                Term {
                    coefficient: 1,
                    exponents: vec![1, 2] // xy^2
                },
                Term {
                    coefficient: modulus - 1,
                    exponents: vec![0, 1] // -y
                }
            ]);
            // basis and divisor
            println!("Basis Polynomials:");
            for poly in &basis {
                println!("{:?}", poly);
            }
            println!("Divisor Polynomial: {:?}", divisor);
            let reduced = divisor.reduce(&basis, modulus);
            // print reduced
            println!("Reduced polynomial: {:?}", reduced);

        }
        else if test == 5 {
            let mut seen_hashes: HashSet<Polynomial> = HashSet::new();
            let mut hasher1 = std::collections::hash_map::DefaultHasher::new();
            p1.hash(&mut hasher1);
            let hash1 = hasher1.finish();
            seen_hashes.insert(p1.clone());

            let mut hasher2 = std::collections::hash_map::DefaultHasher::new();
            p2.hash(&mut hasher2);
            let hash2 = hasher2.finish();
            seen_hashes.insert(p2.clone());

            println!("Hash of Polynomial 1: {}", hash1);
            println!("Hash of Polynomial 2: {}", hash2);
            //print seen_hashes
            println!("Seen Hashes: {:?}", seen_hashes);

            // Create another polynomial identical to p1 to test equality and hashing
            let p1_clone = Polynomial::new(vec![
                Term {
                    coefficient: 1,
                    exponents: vec![2, 0],
                },
                Term {
                    coefficient: modulus - 1,
                    exponents: vec![0, 1],
                },
            ]);
            let mut hasher3 = std::collections::hash_map::DefaultHasher::new();
            p1_clone.hash(&mut hasher3);
            let hash3 = hasher3.finish();
            seen_hashes.insert(p1_clone.clone());

            println!("Hash of Cloned Polynomial 1: {}", hash3);
            println!("Seen Hashes after inserting clone: {:?}", seen_hashes);

            assert_eq!(p1, p1_clone);
            assert_eq!(hash1, hash3);
        }
        else {
            // x + y + z
            let q1 = Polynomial::new(vec![
                Term {
                    coefficient: 1,
                    exponents: vec![1, 0, 0], // x
                },
                Term {
                    coefficient: 1,
                    exponents: vec![0, 1, 0], // y
                },
                Term {
                    coefficient: 1,
                    exponents: vec![0, 0, 1], // z
                },
            ]);
            // xy + xz + yz
            let q2 = Polynomial::new(vec![
                Term {
                    coefficient: 1,
                    exponents: vec![1, 1, 0], // xy
                },
                Term {
                    coefficient: 1,
                    exponents: vec![1, 0, 1], // xz
                },
                Term {
                    coefficient: 1,
                    exponents: vec![0, 1, 1], // yz 
                },
            ]);
            // xyz - 1
            let q3 = Polynomial::new(vec![
                Term {
                    coefficient: 1,
                    exponents: vec![1, 1, 1], // xyz
                },
                Term {
                    coefficient: modulus - 1,
                    exponents: vec![0, 0, 0], // -1
                },
            ]); 
            let p1clone = p1.clone();
            let basis = naive_grobner_basis(vec![p1,p2,p3,p4,p5], modulus);
            // copy basis
            let mut copied_basis = basis.clone();
            println!("Final Grobner Basis:");
            for poly in basis {
                println!("{:?}\n", poly);
            }
            

        /*     // x0 + x1 + x2 + x3
            let test_poly = Polynomial::new(vec![
                Term { coefficient: 1, exponents: vec![1, 0, 0, 0] },
                Term { coefficient: 1, exponents: vec![0, 1, 0, 0] },
                Term { coefficient: 1, exponents: vec![0, 0, 1, 0] },
                Term { coefficient: 1, exponents: vec![0, 0, 0, 1] },
            
            ]);

            // x1^2 + 2x1*x3 +x4^2
            let test_poly_2 = Polynomial::new(vec![
                Term {
                    coefficient: 1,
                    exponents: vec![0, 2, 0, 0],
                },
                Term {
                    coefficient: 2,
                    exponents: vec![0, 1, 0, 1],
                },

                Term {
                    coefficient: 1,
                    exponents: vec![0, 0, 0, 2],
                },
            
            ]);

            // x1*x2 - x1*x3 + x2^2*x3^4 + x2*x3 - 2*x3^2
            let test_poly_3 = Polynomial::new(vec![
            Term {
                coefficient: 1,
                exponents: vec![0, 1, 1, 0], 
            },
            Term {
                coefficient: modulus - 1,
                exponents: vec![0, 1, 0, 1], 
            },
            Term {
                coefficient: 1,
                exponents: vec![0, 0, 2, 4], 
            },
            Term {
                coefficient: 1,
                exponents: vec![0, 0, 1, 1]
            },
            Term {
                coefficient: modulus - 2,
                exponents: vec![0, 0, 0, 2]
            }]);

            // x1*x3^4 - x1 +x3^5 -x3
            let test_poly_4 = Polynomial::new(vec![
                Term{ coefficient: 1, exponents: vec![0, 1, 0, 4] },
                Term{ coefficient: modulus - 1, exponents: vec![0, 1, 0, 0] },
                Term{ coefficient: 1, exponents: vec![0, 0, 0, 5] },
                Term{ coefficient: modulus - 1, exponents: vec![0, 0, 0, 1] }
            ]);
            
            // x2^3*x3^2 + x2^2*x3^3 - x2 -x3
            let test_poly_5 = Polynomial::new(vec![
                Term{ coefficient: 1, exponents: vec![0, 0, 3, 2] },
                Term{ coefficient: 1, exponents: vec![0, 0, 2, 3] },
                Term{ coefficient: modulus - 1, exponents: vec![0, 0, 1, 0] },
                Term{ coefficient: modulus - 1, exponents: vec![0, 0, 0, 1] }
            ]);

            //x2^2*x3^6 - x2^2*x3^2 -x3^4 + 1
            let test_poly_6 = Polynomial::new(vec![
                Term{ coefficient: 1, exponents: vec![0, 0, 2, 6] },
                Term{ coefficient: modulus - 1, exponents: vec![0, 0, 2, 2] },
                Term{ coefficient: modulus - 1, exponents: vec![0, 0, 0, 4] },
                Term{ coefficient: 1, exponents: vec![0, 0, 0, 0] }
            ]);

           
            // x0 + x1 + x2
            let testq1 = Polynomial::new(vec![
                Term{ coefficient: 1, exponents: vec![1, 0, 0, 0] },
                Term{ coefficient: 1, exponents: vec![0, 1, 0, 0] },
                Term{ coefficient: 1, exponents: vec![0, 0, 1, 0] }
            ]);

            // x0^2 + x0x1 + x1^2
            let testq2 = Polynomial::new(vec![
                Term{ coefficient: 1, exponents: vec![2, 0, 0, 0] },
                Term{ coefficient: 1, exponents: vec![1, 1, 0, 0] },
                Term{ coefficient: 1, exponents: vec![0, 2, 0, 0] }
            ]);

            //x0^3 + 32002
            let testq3 = Polynomial::new(vec![
                Term{ coefficient: 1, exponents: vec![3, 0, 0, 0] },
                Term{ coefficient: modulus - 1, exponents: vec![0, 0, 0, 0] }
            ]);


        
            println!("eh");
            let test_basis = vec![test_poly, test_poly_2, test_poly_3, test_poly_4, test_poly_5, test_poly_6];
            let test_basis_clone = test_basis.clone();
            let test = copied_basis[1].clone();
            let is_equivalent = are_bases_equivalent(copied_basis, test_basis, modulus);
            println!("Are the computed basis and test basis equivalent? {}", is_equivalent);


            // print test_basis_clone
            println!("Test Basis Polynomials:");
            for poly in &test_basis_clone {
                println!("{:?}\n", poly);
            }
 */

        }
    }
}