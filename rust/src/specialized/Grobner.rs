use std::collections::HashSet;

use std::hash::{Hash, Hasher};

#[derive(Clone, Debug, PartialEq)]
pub struct Term {
    pub coefficient: f64,
    pub exponents: Vec<usize>, // Exponents for each variable
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
        // Sort terms lexicographically by their exponents in descending order (leading term first)
        terms.sort_by(|a, b| b.exponents.cmp(&a.exponents));
        terms.retain(|t| t.coefficient != 0.0); // Remove zero coefficient terms
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
                            println!("Leading Term: {:?}, Divisor Leading Term: {:?}", leading_term, divisor_leading_term);

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

        result
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
    for i in 0..2 { // This is *supposed* to go until no new polynomials are added, but for now just do 3 iterations
        let mut new_basis = Vec::new();
        for i in 0..basis.len() {
            for j in i + 1..basis.len() {
                let s_poly = Polynomial::s_polynomial(&basis[i], &basis[j]);
                let reduced = s_poly.reduce(&basis);
                //print basis[i], basis[j], s_poly, reduced
                println!("Basis 1: {:?}", basis[i]);
                println!("Basis 2: {:?}", basis[j]);
                println!("S-Polynomial: {:?}", s_poly);
                println!("Reduced: {:?}", reduced);
                if !reduced.terms.is_empty() && !basis_set.contains(&reduced) {
                    println!("Adding new polynomial to basis.");
                    new_basis.push(reduced.clone());
                    basis_set.insert(reduced);
                }
                else {
                    println!("Reduced polynomial is zero or already in basis, skipping.");
                }
            }
        }

        if new_basis.is_empty() {
            break;
        }

        basis.extend(new_basis);

        //print basis with new lines separating each polynomial
        println!("New basis polynomials:");
        for poly in &basis {
            println!("{:?}", poly);
        }
        println!("End of iteration {}\n", i);
    }

    basis
}

fn main() {
    // 1 for s_polynomial, 2 for add, 3 for subtract, 4 for reduce, 5 for testing hashes, else grobner basis
    let test = 0;
    // x^2*y + y^2*z + x*z^2
    let p1 = Polynomial::new(vec![
        Term {
            coefficient: 1.0,
            exponents: vec![2, 1, 0],
        },
        Term {
            coefficient: 1.0,
            exponents: vec![0, 2, 1],
        },
        Term {
            coefficient: 1.0,
            exponents: vec![1, 0, 2],
        },
    ]);

    // xyz-1
    let p2 = Polynomial::new(vec![
        Term {
            coefficient: 1.0,
            exponents: vec![1, 1, 1],
        },
        Term {
            coefficient: -1.0,
            exponents: vec![0, 0, 0],
        },
    ]);

    // x+y+z
    let p3 = Polynomial::new(vec![
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
        // x^2 - y
        let p1 = Polynomial::new(vec![
            Term {
                coefficient: 1.0,
                exponents: vec![2, 0], // x^2
            },
            Term {
                coefficient: -1.0,
                exponents: vec![0, 1], // -y
            },
        ]);
        // xy - 1
        let p2 = Polynomial::new(vec![
            Term {
                coefficient: 1.0,
                exponents: vec![1, 1], // xy
            },
            Term {
                coefficient: -1.0,
                exponents: vec![0, 0], // -1
            },
        ]);
        // x - y^2
        let p3 = Polynomial::new(vec![
            Term {
                coefficient: 1.0,
                exponents: vec![1, 0], // x
            },
            Term {
                coefficient: -1.0,
                exponents: vec![0, 2], // -y^2
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
        println!("Final Grobner Basis:");
        for poly in basis {
            println!("{:?}", poly);
        }
    }
}