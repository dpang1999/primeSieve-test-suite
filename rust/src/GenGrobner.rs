use core::fmt;
use crate::generic::i_field::IField;
use crate::generic::i_exponent::IExponent;
pub mod generic;
use crate::generic::double_field::DoubleField;
use crate::generic::int_mod_p::IntModP;
use crate::generic::vec_exponent::VecExponent;
use crate::generic::bit_packed_exponent::BitPackedExponent;
pub mod helpers;
use crate::helpers::lcg::lcg;

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

#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub struct Polynomial<C, E>
where
    C: IField + Clone + Hash + Eq + fmt::Debug,
    E: IExponent + Clone + Hash + Eq + fmt::Debug,
{
    pub terms: Vec<Term<C, E>>, // Generic terms
}

impl<C, E> Polynomial<C, E>
where
    C: IField + Clone + Hash + Eq + fmt::Debug,
    E: IExponent + Clone + Hash + Eq + fmt::Debug,
{
    pub fn new(mut terms: Vec<Term<C, E>>) -> Self {
        terms.sort_by(|a, b| b.compare(a));
        terms.retain(|t| !t.coefficient.is_zero());
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
                    coefficient: term.coefficient.clone().s(&term.coefficient.zero()),
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
    C: IField + Clone + Hash + Eq + fmt::Debug,
    E: IExponent + Clone + Hash + Eq + fmt::Debug,
{
    let mut basis = polynomials.clone();
    let mut basis_set: HashSet<Polynomial<C, E>> = HashSet::new();

    for poly in &basis {
        basis_set.insert(poly.clone());
    }

    loop {
        let mut added = false;
        let basis_len = basis.len();

        for i in 0..basis_len {
            for j in i + 1..basis_len {
                let s_poly = Polynomial::s_polynomial(&basis[i], &basis[j]);
                let reduced = s_poly.reduce(&basis);

                if !reduced.terms.is_empty() && !basis_set.contains(&reduced) {
                    basis_set.insert(reduced.clone());
                    basis.push(reduced);
                    println!("Added new polynomial to basis, total size now: {}", basis.len());
                    added = true;
                }
            }
        }

        if !added {
            break;
        }
    }
    basis
}
fn main() {
    println!("This is a generic Grobner basis computation module.");
    TERM_ORDER.set(TermOrder::Lex).expect("TERM_ORDER already initialized");
    let p1 = Polynomial::new(vec![
        Term::from_exponents(DoubleField::new(1.0), VecExponent::new(vec![2, 0])),
        Term::from_exponents(DoubleField::new(-1.0), VecExponent::new(vec![0, 1])),
    ]);
    let p2 = Polynomial::new(vec![
        Term::from_exponents(DoubleField::new(1.0), VecExponent::new(vec![1, 1])),
        Term::from_exponents(DoubleField::new(-1.0), VecExponent::new(vec![0, 0])),
    ]);
    let basis = naive_grobner_basis(vec![p1, p2]);
    println!("Computed Grobner basis with {} polynomials.", basis.len());
    // print basis
    for (i, poly) in basis.iter().enumerate() {
        println!("Polynomial {}: {:?}", i + 1, poly);
    }
}