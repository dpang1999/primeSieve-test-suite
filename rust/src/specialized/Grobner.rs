use std::collections::HashSet;

use std::hash::{Hash, Hasher};

use std::sync::OnceLock;

static TERM_ORDER: OnceLock<TermOrder> = OnceLock::new();

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
    for poly in &basis {
        println!("{:?}", poly);
        basis_set.insert(poly.clone());
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
    TERM_ORDER.set(TermOrder::GrLex).expect("TERM_ORDER already initialized");

    // x^2y + y^2z + z^2x
    let p1 = Polynomial::new(vec![
        Term {
            coefficient: 1.0,
            exponents: vec![2, 1, 0], // x^2 * y
        },
        Term {
            coefficient: 1.0,
            exponents: vec![0, 2, 1], // y^2 * z
        },
        Term {
            coefficient: 1.0,
            exponents: vec![1, 0, 2], // z^2 * x
        },
        

    ]);
    // x*y*z -1
    let p2 = Polynomial::new(vec![
        Term {
            coefficient: 1.0,
            exponents: vec![1, 1, 1], // x*y*z
        },
        Term {
            coefficient: -1.0,
            exponents: vec![0, 0, 0], // -1
        }
    ]);
    // x+y+z
    let p3 = Polynomial::new(vec![
        Term {
            coefficient: 1.0,
            exponents: vec![1, 0, 0], // x
        },
        Term {
            coefficient: 1.0,
            exponents: vec![0, 1, 0], // y
        },
        Term {
            coefficient: 1.0,
            exponents: vec![0, 0, 1], // z
        },
    ]);

    // x + y + z
    let p4 = Polynomial::new(vec![
        Term {
            coefficient: 1.0,
            exponents: vec![1, 0, 0],
        },
        Term {
            coefficient: 1.0,
            exponents: vec![0, 1, 0],
        },
        Term {
            coefficient: 1.0,
            exponents: vec![0, 0, 1],
        },
    ]);

    println!("Begin the experiment");


    if test == 1 {
        let s_poly = Polynomial::s_polynomial(&p1, &p2);

        let expected = Polynomial::new(vec![
            Term {
                coefficient: -1.0,
                exponents: vec![0, 2], // Term: -y^2
            },
            Term {
                coefficient: 1.0,
                exponents: vec![1, 0], // Term: x
            },
        ]);

        assert_eq!(s_poly, expected);
    } else if test == 2 {
        let sum = p1.add(&p2);
        println!("Sum: {:?}", sum);

        let expected = Polynomial::new(vec![
            Term {
                coefficient: 1.0,
                exponents: vec![2, 0],
            },
            Term {
                coefficient: 1.0,
                exponents: vec![1, 1],
            },
            Term {
                coefficient: -1.0,
                exponents: vec![0, 1],
            },
            Term {
                coefficient: -1.0,
                exponents: vec![0, 0],
            },
        ]);

        assert_eq!(sum, expected);
    } else if test == 3 {
        let difference = p1.subtract(&p2);
        println!("Difference: {:?}", difference);

        let expected = Polynomial::new(vec![
            Term {
                coefficient: 1.0,
                exponents: vec![2, 0],
            },
            Term {
                coefficient: -1.0,
                exponents: vec![1, 1],
            },
            Term {
                coefficient: -1.0,
                exponents: vec![0, 1],
            },
            Term {
                coefficient: 1.0,
                exponents: vec![0, 0],
            },
        ]);

        assert_eq!(difference, expected);
    }
    else if test == 4 {
        // x^3 + y^3 + z^3
        let p1 = Polynomial::new(vec![
            Term {
                coefficient: 1.0,
                exponents: vec![3, 0, 0], // x^3
            },
            Term {
                coefficient: 1.0,
                exponents: vec![0, 3, 0], // y^3
            },
            Term {
                coefficient: 1.0,
                exponents: vec![0, 0, 3], // z^3
            },
        ]);
        // xy + yz + xz
        let p2 = Polynomial::new(vec![
            Term {
                coefficient: 1.0,
                exponents: vec![1, 1, 0], // xy
            },
            Term {
                coefficient: 1.0,
                exponents: vec![0, 1, 1], // yz
            },
            Term {
                coefficient: 1.0,
                exponents: vec![1, 0, 1], // xz
            },
        ]);
        // x+y+z
        let p3 = Polynomial::new(vec![
            Term {
                coefficient: 1.0,
                exponents: vec![1, 0, 0], // x
            },
            Term {
                coefficient: 1.0,
                exponents: vec![0, 1, 0], // y
            },
            Term {
                coefficient: 1.0,
                exponents: vec![0, 0, 1], // z
            },
        ]);
        let basis = vec![p1.clone(), p2.clone(), p3.clone()];

        // xy^2 - y
        let divisor = Polynomial::new(vec![
            Term {
                coefficient: 1.0,
                exponents: vec![1, 2] // xy^2
            },
            Term {
                coefficient: -1.0,
                exponents: vec![0, 1] // -y
            }
        ]);
        // basis and divisor
        println!("Basis Polynomials:");
        for poly in &basis {
            println!("{:?}", poly);
        }
        println!("Divisor Polynomial: {:?}", divisor);
        let reduced = divisor.reduce(&basis);
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
                coefficient: 1.0,
                exponents: vec![2, 0],
            },
            Term {
                coefficient: -1.0,
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

        let basis = naive_grobner_basis(vec![p1, p2, p3]);
        // copy basis
        let mut copied_basis = basis.clone();
        println!("Final Grobner Basis:");
        for poly in basis {
            println!("{:?}", poly);
        }
        

        //  x^2 + y^2 -1
        let test_poly = Polynomial::new(vec![
            Term {
                coefficient: 1.0,
                exponents: vec![2, 0, 0],
            },
            Term {
                coefficient: 1.0,
                exponents: vec![0, 2, 0],
            },
            Term {
                coefficient: -1.0,
                exponents: vec![0, 0, 0],
            },
        ]);

        // x*y -z
        let test_poly_2 = Polynomial::new(vec![
            Term {
                coefficient: 1.0,
                exponents: vec![1, 1, 0],
            },
            Term {
                coefficient: -1.0,
                exponents: vec![0, 0, 1],
            },
           
        ]);

        // x*z +y^3 - y
        let test_poly_3 = Polynomial::new(vec![
            Term {
                coefficient: 1.0,
                exponents: vec![1, 0, 1],
            },
            Term {
                coefficient: 1.0,
                exponents: vec![0, 3, 0],
            },
            Term {
                coefficient: -1.0,
                exponents: vec![0, 1, 0],
            },
        ]);

        // y^4 - y^2 + z^2
        let test_poly_4 = Polynomial::new(vec![
            Term {
                coefficient: 1.0,
                exponents: vec![0, 4, 0],
            },
            Term {
                coefficient: -1.0,
                exponents: vec![0, 2, 0],
            },
            Term {
                coefficient: 1.0,
                exponents: vec![0, 0, 2],
            },
        ]);

        let test_basis = vec![test_poly, test_poly_2, test_poly_3];

        // print test_basis
        println!("Test Basis Polynomials:");
        for poly in &test_basis {
            println!("{:?}", poly);
        }
        let is_equivalent = are_bases_equivalent(copied_basis, test_basis);
        println!("Are the computed basis and test basis equivalent? {}", is_equivalent);
    }
}