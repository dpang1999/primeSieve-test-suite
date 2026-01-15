use core::fmt;
use std::env;
use crate::generic::i_field::IField;
use crate::generic::i_exponent::IExponent;
pub mod generic;
use crate::generic::double_field::DoubleField;
use crate::generic::int_mod_p::IntModP;
use crate::generic::single_field::SingleField;
use crate::generic::vec_exponent::VecExponent;
use crate::generic::bit_packed_exponent::BitPackedExponent;
pub mod helpers;
use crate::helpers::lcg::Lcg;

use std::collections::HashSet;
use std::hash::{Hash, Hasher};

use std::sync::OnceLock;

static TERM_ORDER: OnceLock<TermOrder> = OnceLock::new();

#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub struct Term<C, E> 
where 
    C: IField,
    E: IExponent
{
    pub coefficient: C, // Generic coefficient
    pub exponents: E,    // Generic exponents
}

#[derive(Copy, Clone, Debug, PartialEq)]
pub enum TermOrder {
    Lex,
    GrLex,
    RevLex
}


impl<C,E> Term<C, E>
where 
    C: IField,
    E: IExponent
{
    pub fn from_exponents(coefficient: C, exponents: E) -> Self {
        Term { coefficient, exponents }
    }

    pub fn compare(&self, other: &Term<C, E>) -> std::cmp::Ordering {
        match TERM_ORDER.get().unwrap() {
            TermOrder::Lex => self.exponents.lex_compare(&other.exponents),
            TermOrder::GrLex => {
                let degree_self = self.exponents.degree();
                let degree_other = other.exponents.degree();
                if degree_self != degree_other {
                    degree_self.cmp(&degree_other)
                } else {
                    self.exponents.lex_compare(&other.exponents)
                }
            }
            TermOrder::RevLex => {
                let degree_self = self.exponents.degree();
                let degree_other = other.exponents.degree();
                if degree_self != degree_other {
                    degree_self.cmp(&degree_other)
                } else {
                    other.exponents.lex_compare(&self.exponents)
                }
            }
        }
    }

    pub fn can_reduce(&self, divisor_leading: &Term<C, E>) -> bool {
        self.exponents.can_reduce(&divisor_leading.exponents)
    }

    pub fn lcm(&self, other: &Term<C, E>) -> E {
        self.exponents.lcm(&other.exponents)
    }
}

impl<C, E> fmt::Display for Term<C, E>
where 
    C: IField + fmt::Display,
    E: IExponent + fmt::Display
{
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}*{}", self.coefficient, self.exponents)
    }
}

#[derive(Clone, PartialEq, Eq, Hash)]
pub struct Polynomial<C, E>
where
    C: IField + Clone + Hash + Eq + fmt::Display,
    E: IExponent + Clone + Hash + Eq + fmt::Display,
{
    pub terms: Vec<Term<C, E>>, // Generic terms
}

impl<C, E> fmt::Display for Polynomial<C, E>
where
    C: IField + Clone + Hash + Eq + fmt::Display,
    E: IExponent + Clone + Hash + Eq + fmt::Display,
{
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{{ ")?;
        for (i, term) in self.terms.iter().enumerate() {
            if i > 0 {
                write!(f, " + ")?;
            }
            write!(f, "{}*{}", term.coefficient, term.exponents)?;
        }
        write!(f, " }}")
    }
}

impl<C, E> Polynomial<C, E>
where
    C: IField + Clone + Hash + Eq + fmt::Display,
    E: IExponent + Clone + Hash + Eq + fmt::Display,
{
    pub fn new(mut terms: Vec<Term<C, E>>) -> Self {
        terms.sort_by(|a, b| b.compare(a));
        //terms.retain(|t| !t.coefficient.is_zero());
        terms.retain(|t| (t.coefficient.coerce_to_f64() - 0.0).abs() > 1e-2);
        Polynomial { terms }
    }

    pub fn add(&self, other: &Polynomial<C, E>) -> Polynomial<C, E> {
        let mut result = self.terms.clone();
        for term in &other.terms {
            let mut found = false;
            for res_term in &mut result {
                if res_term.exponents == term.exponents {
                    res_term.coefficient = res_term.coefficient.a(&term.coefficient);
                    found = true;
                    break;
                }
            }
            if !found {
                result.push(term.clone());
            }
        }
        Polynomial::new(result)
    }

    pub fn subtract(&self, other: &Polynomial<C, E>) -> Polynomial<C, E> {
        let mut result = self.terms.clone();
        for term in &other.terms {
            let mut found = false;
            for res_term in &mut result {
                if res_term.exponents == term.exponents {
                    res_term.coefficient = res_term.coefficient.s(&term.coefficient);
                    found = true;
                    break;
                }
            }
            if !found {
                let neg_term = Term {
                    coefficient: term.coefficient.zero().s(&term.coefficient.clone()),
                    exponents: term.exponents.clone(),
                };
                result.push(neg_term);
            }
        }
        Polynomial::new(result)
    }

    pub fn reduce(&self, divisors: &[Polynomial<C, E>]) -> Polynomial<C, E> {
        let mut result = self.clone();

        loop {
            let mut reduced = false;

            for divisor in divisors {
                if let Some(leading_term) = result.terms.first() {
                    if let Some(divisor_leading_term) = divisor.terms.first() {
                        if leading_term.can_reduce(divisor_leading_term) {
                            // debug
                            //println!("Reducing term: {} by divisor leading term: {}", leading_term, divisor_leading_term);
                            let coefficient = leading_term
                                .coefficient
                                .d(&divisor_leading_term.coefficient);
                            let exponents = leading_term.exponents.sub(&divisor_leading_term.exponents);

                            let reduction_term = Term {
                                coefficient,
                                exponents,
                            };

                            let scaled_divisor = divisor.multiply_by_term(&reduction_term);
                            result = result.subtract(&scaled_divisor);
                            //println!("After reduction, polynomial is: {}", result);
                            reduced = true;
                            break;
                        }
                    }
                }
            }

            if !reduced {
                break;
            }
        }

        Polynomial::new(result.terms)
    }

    pub fn multiply_by_term(&self, term: &Term<C, E>) -> Polynomial<C, E> {
        let terms = self
            .terms
            .iter()
            .map(|t| Term {
                coefficient: t.coefficient.m(&term.coefficient),
                exponents: t.exponents.add(&term.exponents),
            })
            .collect();

        Polynomial::new(terms)
    }

    pub fn s_polynomial(p1: &Polynomial<C, E>, p2: &Polynomial<C, E>) -> Polynomial<C, E> {
        let leading_term_p1 = &p1.terms[0];
        let leading_term_p2 = &p2.terms[0];

        let lcm_exponents = leading_term_p1.lcm(&leading_term_p2);

        let scale_factor_p1 = lcm_exponents.sub(&leading_term_p1.exponents);
        let scale_factor_p2 = lcm_exponents.sub(&leading_term_p2.exponents);



        let scaled_p1 = p1.multiply_by_term(&Term{
            coefficient: leading_term_p1.coefficient.one(),
            exponents: scale_factor_p1
        });

        let scaled_p2 = p2.multiply_by_term(&Term{
            coefficient: leading_term_p2.coefficient.one(),
            exponents: scale_factor_p2
        });
        scaled_p1.subtract(&scaled_p2)
    }
}

pub fn naive_grobner_basis<C, E>(polynomials: Vec<Polynomial<C, E>>) -> Vec<Polynomial<C, E>>
where
    C: IField + Clone + Hash + Eq + fmt::Display,
    E: IExponent + Clone + Hash + Eq + fmt::Display,
{
    let mut basis = polynomials.clone();
    let mut basis_set: HashSet<Polynomial<C, E>> = HashSet::new();

    for poly in &basis {
        //println!("Initial basis polynomial: {}", poly);
        basis_set.insert(poly.clone());
    }

    let mut polyAdded = 0;
    loop {
        polyAdded += 1;
        println!("Grobner basis iteration {}", polyAdded);
        let mut added = false;
        let basis_len = basis.len();

        for i in 0..basis_len {
            for j in i + 1..basis_len {
                let s_poly = Polynomial::s_polynomial(&basis[i], &basis[j]);
                let reduced = s_poly.reduce(&basis);
                //print basis i and basis j

                if !reduced.terms.is_empty() && !basis_set.contains(&reduced) {
                    //println!("Reducing S-Polynomial of basis polynomials {} and {}", &basis[i], &basis[j]);
                    //println!("Reduced S-Polynomial of basis[{}] and basis[{}]: {}", i, j, &reduced);
                    basis_set.insert(reduced.clone());
                    basis.push(reduced);
                    //println!("Added new polynomial to basis, total size now: {}", basis.len());
                    added = true;
                }
            }
        }

        if !added {
            break;
        }
    }

    // print basis before reduction
    /*println!("Basis before reduction:");
    for poly in &basis {
        println!("{}", poly);
    }*/

    let mut reduced_basis = Vec::new();
    for poly in &basis {
        let mut basis_excluding_self = basis.clone();
        basis_excluding_self.retain(|p| p != poly);
        let reduced = poly.reduce(&basis_excluding_self);
        if !reduced.terms.is_empty() && !reduced_basis.contains(&reduced) {
            reduced_basis.push(reduced);
        }
    }
    reduced_basis
}

fn main() {
    // let mode = 0 be for testing
    let mode = 1;
    println!("This is a generic Grobner basis computation module.");
    if mode != 0 {
        let mut rand = Lcg::new(12345, 1345, 65, 17);
        let args: Vec<String> = env::args().collect();

        // arg1 = # of polynomials
        // arg2 = coefficient type (0 = SingleField, 1 = DoubleField, 2 = IntModP)
        // arg3 = exponent type (0 = VecExponent, 1 = BitPackedExponent) 
        // arg4 = term order (0 = Lex, 1 = GrLex, 2 = RevLex)
        let polynomial_count = if args.len() > 1 {
            args[1].parse::<usize>().unwrap_or(3)
        } else {
            3
        };
        let coeff_type = if args.len() > 2 {
            args[2].parse::<usize>().unwrap_or(0)
        } else {
            0
        };
        let exp_type = if args.len() > 3 {
            args[3].parse::<usize>().unwrap_or(1)
        } else {
            1
        };
        let term_order = if args.len() > 4 {
            args[4].parse::<usize>().unwrap_or(0)
        } else {
            0
        };
        match term_order {
            0 => { TERM_ORDER.set(TermOrder::Lex).expect("TERM_ORDER already initialized"); },
            1 => { TERM_ORDER.set(TermOrder::GrLex).expect("TERM_ORDER already initialized"); },
            2 => { TERM_ORDER.set(TermOrder::RevLex).expect("TERM_ORDER already initialized"); },
            _ => { TERM_ORDER.set(TermOrder::Lex).expect("TERM_ORDER already initialized"); },
        }

        let generated = generate_polynomials(polynomial_count, coeff_type, exp_type, &mut rand);

        match generated {
            GeneratedPolynomials::SingleFieldVecExponent(polys) => {
                let basis = naive_grobner_basis(polys);
                println!("Computed Grobner basis with {} polynomials.", basis.len());
                /*for (i, poly) in basis.iter().enumerate() {
                    println!("Polynomial {}: {}", i + 1, poly);
                }*/
            }
            GeneratedPolynomials::SingleFieldBitPackedExponent(polys) => {
                let basis = naive_grobner_basis(polys);
                println!("Computed Grobner basis with {} polynomials.", basis.len());
                /*for (i, poly) in basis.iter().enumerate() {
                    println!("Polynomial {}: {}", i + 1, poly);
                }*/
            }
            GeneratedPolynomials::DoubleFieldVecExponent(polys) => {
                let basis = naive_grobner_basis(polys);
                println!("Computed Grobner basis with {} polynomials.", basis.len());
                /*for (i, poly) in basis.iter().enumerate() {
                    println!("Polynomial {}: {}", i + 1, poly);
                }*/
            }
            GeneratedPolynomials::DoubleFieldBitPackedExponent(polys) => {
                let basis = naive_grobner_basis(polys);
                println!("Computed Grobner basis with {} polynomials.", basis.len());
                /*for (i, poly) in basis.iter().enumerate() {
                    println!("Polynomial {}: {}", i + 1, poly);
                }*/
            }
            GeneratedPolynomials::IntModPVecExponent(polys) => {
                let basis = naive_grobner_basis(polys);
                println!("Computed Grobner basis with {} polynomials.", basis.len());
                /*for (i, poly) in basis.iter().enumerate() {
                    println!("Polynomial {}: {}", i + 1, poly);
                }*/
            }
            GeneratedPolynomials::IntModPBitPackedExponent(polys) => {
                let basis = naive_grobner_basis(polys);
                println!("Computed Grobner basis with {} polynomials.", basis.len());
                /*for (i, poly) in basis.iter().enumerate() {
                    println!("Polynomial {}: {}", i + 1, poly);
                }*/
            }
            GeneratedPolynomials::None => {
                println!("No polynomials generated.");
            }
        }


    }
    else if mode == 0 {
        TERM_ORDER.set(TermOrder::Lex).expect("TERM_ORDER already initialized");
        //x^2*y + y^2*z + z^2*x
        let p1 = Polynomial::new(vec![
            Term::from_exponents(SingleField::new(1.0), BitPackedExponent::from_vec([2,1,0,0,0,0])),
            Term::from_exponents(SingleField::new(1.0), BitPackedExponent::from_vec([0,2,1,0,0,0])),
            Term::from_exponents(SingleField::new(1.0), BitPackedExponent::from_vec([1,0,2,0,0,0])),
        ]);
        //xyz - 1
        let p2 = Polynomial::new(vec![
            Term::from_exponents(SingleField::new(1.0), BitPackedExponent::from_vec([1, 1,1,0,0,0])),
            Term::from_exponents(SingleField::new(-1.0), BitPackedExponent::from_vec([0, 0,0,0,0,0])),
        ]);
        // x+ y + z
        let p3 = Polynomial::new(vec![
            Term::from_exponents(SingleField::new(1.0), BitPackedExponent::from_vec([1,0,0,0,0,0])),
            Term::from_exponents(SingleField::new(1.0), BitPackedExponent::from_vec([0,1,0,0,0,0])),
            Term::from_exponents(SingleField::new(1.0), BitPackedExponent::from_vec([0,0,1,0,0,0])),
        ]);
        let initial_basis = vec![p1.clone(), p2.clone(), p3.clone()];
        for(i, poly) in initial_basis.iter().enumerate() {
            println!("Initial Polynomial {}: {}", i + 1, poly);
        }
        let basis = naive_grobner_basis(vec![p1, p2, p3]);
        println!("Computed Grobner basis with {} polynomials.", basis.len());
        // print basis
        for (i, poly) in basis.iter().enumerate() {
            println!("Polynomial {}: {}", i + 1, poly);
        }
    }
}

// Helper function to generate a vector of polynomials for all type combinations
// Returns a tuple of Option<Vec<Polynomial<...>>> for each type combination
// conditions: always 3 terms per poly, 3 variables, exponents in [0,3]
// coefficient type (0 = SingleField, 1 = DoubleField, 2 = IntModP)
// exponent type (0 = VecExponent, 1 = BitPackedExponent) 
pub fn generate_polynomials(
    num_polys: usize,
    coeff_type: usize,
    exp_type: usize,
    rand: &mut Lcg,
) -> GeneratedPolynomials {
    let max_exp_value = 4; // Exponents will be in the range [0, max_exp_value-1]
    match (coeff_type, exp_type) {
        (0, 0) => GeneratedPolynomials::SingleFieldVecExponent(
            (0..num_polys)
                .map(|_| {
                    let terms = (0..3).map(|_| {
                        let coeff = SingleField::new(rand.next_double() as f32);
                        let exps = vec![(rand.next_int() % max_exp_value) as u32, (rand.next_int() % max_exp_value) as u32, (rand.next_int() % max_exp_value) as u32];
                        Term::from_exponents(coeff, VecExponent::new(exps))
                    }).collect();
                    Polynomial::new(terms)
                })
                .collect()
        ),
        (0, 1) => GeneratedPolynomials::SingleFieldBitPackedExponent(
            (0..num_polys)
                .map(|_| {
                    let terms = (0..3).map(|_| {
                        let coeff = SingleField::new(rand.next_double() as f32);
                        let mut arr = [0u8; 6];
                        for i in 0..3 { arr[i] = (rand.next_int() % max_exp_value) as u8; }
                        Term::from_exponents(coeff, BitPackedExponent::from_vec(arr))
                    }).collect();
                    Polynomial::new(terms)
                })
                .collect()
        ),
        (1, 0) => GeneratedPolynomials::DoubleFieldVecExponent(
            (0..num_polys)
                .map(|_| {
                    let terms = (0..3).map(|_| {
                        let coeff = DoubleField::new(rand.next_double());
                        let exps = vec![(rand.next_int() % max_exp_value) as u32, (rand.next_int() % max_exp_value) as u32, (rand.next_int() % max_exp_value) as u32];
                        Term::from_exponents(coeff, VecExponent::new(exps))
                    }).collect();
                    Polynomial::new(terms)
                })
                .collect()
        ),
        (1, 1) => GeneratedPolynomials::DoubleFieldBitPackedExponent(
            (0..num_polys)
                .map(|_| {
                    let terms = (0..3).map(|_| {
                        let coeff = DoubleField::new(rand.next_double());
                        let mut arr = [0u8; 6];
                        for i in 0..3 { arr[i] = (rand.next_int() % max_exp_value) as u8; }
                        Term::from_exponents(coeff, BitPackedExponent::from_vec(arr))
                    }).collect();
                    Polynomial::new(terms)
                })
                .collect()
        ),
        (2, 0) => GeneratedPolynomials::IntModPVecExponent(
            (0..num_polys)
                .map(|_| {
                    let terms = (0..3).map(|_| {
                        let coeff = IntModP::new((rand.next_int() % 7) as u128, 7);
                        let exps = vec![(rand.next_int() % max_exp_value) as u32, (rand.next_int() % max_exp_value) as u32, (rand.next_int() % max_exp_value) as u32];
                        Term::from_exponents(coeff, VecExponent::new(exps))
                    }).collect();
                    Polynomial::new(terms)
                })
                .collect()
        ),
        (2, 1) => GeneratedPolynomials::IntModPBitPackedExponent(
            (0..num_polys)
                .map(|_| {
                    let terms = (0..3).map(|_| {
                        let coeff = IntModP::new((rand.next_int() % 7) as u128, 7);
                        let mut arr = [0u8; 6];
                        for i in 0..3 { arr[i] = (rand.next_int() % max_exp_value) as u8; }
                        Term::from_exponents(coeff, BitPackedExponent::from_vec(arr))
                    }).collect();
                    Polynomial::new(terms)
                })
                .collect()
        ),
        _ => GeneratedPolynomials::None,
    }
}

// Enum to represent all possible generated polynomial types
pub enum GeneratedPolynomials {
    SingleFieldVecExponent(Vec<Polynomial<SingleField, VecExponent>>),
    SingleFieldBitPackedExponent(Vec<Polynomial<SingleField, BitPackedExponent>>),
    DoubleFieldVecExponent(Vec<Polynomial<DoubleField, VecExponent>>),
    DoubleFieldBitPackedExponent(Vec<Polynomial<DoubleField, BitPackedExponent>>),
    IntModPVecExponent(Vec<Polynomial<IntModP, VecExponent>>),
    IntModPBitPackedExponent(Vec<Polynomial<IntModP, BitPackedExponent>>),
    None,
}